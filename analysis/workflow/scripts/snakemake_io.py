import os
import re
import collections
from itertools import chain

# Copied from https://github.com/snakemake/snakemake/blob/27b224ed/snakemake/io.py#L828C5-L843
_wildcard_regex = re.compile(
    r"""
    \{
        (?=(   # This lookahead assertion emulates an 'atomic group'
               # which is required for performance
            \s*(?P<name>\w+)                    # wildcard name
            (\s*,\s*
                (?P<constraint>                 # an optional constraint
                    ([^{}]+ | \{\d+(,\d+)?\})*  # allow curly braces to nest one level
                )                               # ...  as in '{w,a{3,5}}'
            )?\s*
        ))\1
    \}
    """,
    re.VERBOSE,
)

def regex(filepattern):
    """
    Copied from https://github.com/snakemake/snakemake/blob/27b224ed/snakemake/io.py#L935C1-L960C22
    """
    f = []
    last = 0
    wildcards = set()
    for match in _wildcard_regex.finditer(filepattern):
        f.append(re.escape(filepattern[last : match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError(
                    "Constraint regex must be defined only in the first "
                    "occurence of the wildcard in a string."
                )
            f.append(f"(?P={wildcard})")
        else:
            wildcards.add(wildcard)
            f.append(
                "(?P<{}>{})".format(
                    wildcard,
                    match.group("constraint") if match.group("constraint") else ".+",
                )
            )
        last = match.end()
    f.append(re.escape(filepattern[last:]))
    f.append("$")  # ensure that the match spans the whole file
    return "".join(f)

def glob_wildcards(pattern, files=None, followlinks=False):
    """
    Glob the values of the wildcards by matching the given pattern to the filesystem.
    Returns a named tuple with a list of values for each wildcard.

    Copied from https://github.com/snakemake/snakemake/blob/27b224ed/snakemake/io.py#L1310C1-L1345C21
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    dirname = (
        os.path.dirname(pattern[: first_wildcard.start()])
        if first_wildcard
        else os.path.dirname(pattern)
    )
    if not dirname:
        dirname = "."

    names = [match.group("name") for match in _wildcard_regex.finditer(pattern)]
    Wildcards = collections.namedtuple("Wildcards", names)
    wildcards = Wildcards(*[list() for name in names])

    pattern = re.compile(regex(pattern))

    if files is None:
        files = (
            os.path.normpath(os.path.join(dirpath, f))
            for dirpath, dirnames, filenames in os.walk(
                dirname, followlinks=followlinks
            )
            for f in chain(filenames, dirnames)
        )

    for f in files:
        match = re.match(pattern, f)
        if match:
            for name, value in match.groupdict().items():
                getattr(wildcards, name).append(value)
    return wildcards
