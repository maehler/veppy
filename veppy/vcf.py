from collections.abc import Iterable
from dataclasses import dataclass
import gzip
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple

from veppy import utils
from veppy.vep import VEP, VEPAnnotation


@dataclass
class Metadata:
    ID: str
    Description: str

    def __post_init__(self):
        self.Description = self.Description.strip("'\"")


@dataclass
class FormatMetadata(Metadata):
    Type: str
    Number: str


@dataclass
class InfoMetadata(Metadata):
    Type: str
    Number: str
    Source: Optional[str] = None
    Version: Optional[str] = None


class Sample:
    def __init__(self, name: str, attr: Dict[str, str]):
        self.name = name
        self._attr_dict = attr

    def __getitem__(self, key: str) -> str:
        return self._attr_dict[key]


class SampleList:
    def __init__(self, samples: Iterable[Sample]):
        self._samples = list(samples)
        sample_names = [s.name for s in samples]
        self._sample_dict = dict(zip(sample_names, self._samples))

    def __getitem__(self, key: str | int) -> Sample:
        if isinstance(key, str):
            return self._sample_dict[key]
        elif isinstance(key, int):
            return self._samples[key]
        else:
            raise TypeError("invalid key, str or int needed")


class Variant:
    def __init__(
        self, variant_string: str, sample_names: List[str], annotation: Optional[VEP]
    ):
        attr = variant_string.split("\t")
        if len(attr) < 10:
            raise ValueError("invalid variant, too few fields")
        self.chrom = attr[0]
        self.pos = int(attr[1])
        self.id = attr[2]
        self.ref = attr[3]
        self.alt = attr[4].split(",")
        self.qual = float(attr[5]) if attr[5] != "." else None
        self.filter = attr[6]
        self.info = self._parse_info(attr[7])
        self.format = self._parse_format(attr[8])
        self.samples = self._parse_samples(attr[9:], sample_names)

        if annotation is not None:
            self._annotation = annotation.variant(self)

    @property
    def annotation(self) -> Optional[VEPAnnotation]:
        return self._annotation

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key.lower())

    def _parse_info(self, info_str: str) -> Dict[str, str]:
        if len(info_str) == 0:
            return dict()
        attr = unquoted_semicolon.split(info_str)
        attr_dict = {}
        for a in attr:
            split_attr = unquoted_equals.split(a)
            if len(split_attr) == 1:
                attr_dict[a] = True
            else:
                attr_dict[split_attr[0]] = split_attr[1]
        return attr_dict

    def _parse_format(self, format_str: str) -> List[str]:
        return format_str.split(":")

    def _parse_samples(self, samples: List[str], names: List[str]) -> SampleList:
        sample_list = []
        for n, s in zip(names, samples):
            sample = Sample(n, dict(zip(self.format, s.split(":"))))
            sample_list.append(sample)
        return SampleList(sample_list)


unquoted_comma = re.compile(r',(?=(?:[^"]*"[^"]*")*[^"]*$)')
unquoted_semicolon = re.compile(r';(?=(?:[^"]*"[^"]*")*[^"]*$)')
unquoted_equals = re.compile(r'=(?=(?:[^"]*"[^"]*")*[^"]*$)')


class VCF:
    def __init__(
        self, path: str | Path, annotation_key: str = "CSQ", annotation_sep: str = "|"
    ) -> None:
        openfun = open
        if utils.is_gzipped(path):
            openfun = gzip.open

        self._annotation_key = annotation_key
        self._annotation_sep = annotation_sep
        self._annotation = None
        self._fp = openfun(path, "rt")
        self._parse()

    @property
    def annotations(self) -> Optional[VEP]:
        return self._annotation

    def _parse_version(self) -> str:
        line = next(self._fp)
        assert line.startswith("##fileformat=")
        return line.strip().split("=")[1]

    def _parse_metadata(
        self,
    ) -> Tuple[List[str], Dict[str, str], Dict[str, str], Optional[VEP]]:
        format = {}
        info = {}
        samples = []
        annotations = None

        for line in self._fp:
            line = line.strip()
            if not line.startswith("##"):
                line = line.split()
                samples = line[9:]
                break

            match unquoted_equals.split(line, maxsplit=1):
                case ["##FORMAT", attr_string]:
                    attr_string = attr_string.strip("<>")
                    attributes = [
                        unquoted_equals.split(a)
                        for a in unquoted_comma.split(attr_string)
                    ]
                    attr_dict = dict((k, v) for k, v in attributes)
                    format[attr_dict["ID"]] = FormatMetadata(**attr_dict)
                case ["##INFO", attr_string]:
                    attr_string = attr_string.strip("<>")
                    attributes = [
                        unquoted_equals.split(a)
                        for a in unquoted_comma.split(attr_string)
                    ]
                    attr_dict = dict((k, v) for k, v in attributes)
                    info[attr_dict["ID"]] = InfoMetadata(**attr_dict)
                    if attr_dict["ID"] == self._annotation_key:
                        annotations = VEP(info[self._annotation_key].Description)
                case other:
                    pass

        return samples, info, format, annotations

    def _parse_variants(self):
        for line in self._fp:
            line = line.strip()
            yield Variant(line, self.samples, self._annotation)

    def _parse(self):
        self.version = self._parse_version()
        self.samples, self.info, self.format, self._annotation = self._parse_metadata()
        self.variants = self._parse_variants()

    def __del__(self):
        self._fp.close()
