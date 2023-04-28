import click
from pathlib import Path
import sys

from veppy.vcf import VCF
from veppy.version import __version__


@click.command()
@click.version_option(version=__version__)
@click.argument(
    "vcf_path",
    metavar="VCF",
    type=click.Path(path_type=Path, exists=True),
)
@click.argument(
    "vep_fields",
    metavar="VEP",
    nargs=-1,
)
@click.option(
    "--list",
    "-l",
    "list_fields",
    is_flag=True,
    help="list annotation fields",
)
@click.option(
    "--fields",
    "-f",
    "vcf_fields",
    metavar="FIELD",
    help="VCF fields to include in output",
    type=click.Choice(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]),
    default=["CHROM", "POS", "REF", "ALT"],
    show_default=True,
    multiple=True,
)
@click.option(
    "--vep-key",
    "-k",
    metavar="KEY",
    help="annotation key in INFO field",
    default="CSQ",
    show_default=True,
)
@click.option(
    "--sep",
    metavar="SEP",
    help="annotation separator",
    default="|",
    show_default=True,
)
def main(vcf_path, vep_fields, list_fields, vcf_fields, vep_key, sep):
    vcf = VCF(vcf_path, annotation_key=vep_key, annotation_sep=sep)

    if vcf.annotations is None:
        print("error: no annotations found", file=sys.stderr)
        sys.exit(1)

    if list_fields:
        print("\n".join(vcf.annotations.variables))
        return

    if len(vep_fields) == 0:
        vep_fields = tuple(vcf.annotations.variables)
    else:
        for f in vep_fields:
            if f not in vcf.annotations.variables:
                print(f"error: annotation field not found: {f}", file=sys.stderr)
                sys.exit(1)

    header = vcf_fields + vep_fields
    print("\t".join(header))

    for v in vcf.variants:
        if v.annotation is None:
            print("warning: no annotations for variant", file=sys.stderr)
            continue
        cols = []
        for f in vcf_fields:
            cols.append(str(v[f]))
        for f in vep_fields:
            if v.annotation[f] is None:
                cols.append("")
            else:
                cols.append(str(v.annotation[f]))
        print("\t".join(cols))
