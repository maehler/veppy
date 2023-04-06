from typing import TYPE_CHECKING, Dict, Optional

if TYPE_CHECKING:
    from veppy.vcf import Variant


VEPAnnotation = Dict[str, Optional[str]]


class VEP:
    def __init__(self, description: str, key: str = "CSQ", sep: str = "|") -> None:
        self._key = key
        self._sep = sep
        self._variables = description.split()[-1].split(self._sep)

    @property
    def variables(self):
        return self._variables

    def variant(self, variant: "Variant") -> Optional[VEPAnnotation]:
        vep_dict = {}
        if self._key not in variant.info:
            return None
        for k, v in zip(self.variables, variant.info[self._key].split(self._sep)):
            vep_dict[k] = v if v != "" else None
        return vep_dict
