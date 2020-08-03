"""XYZ file format readers."""


from pyparsing import (Group, LineEnd, OneOrMore, Regex,
                       Suppress, Word, alphas, nums, restOfLine)


__all__ = ["parser_xyz"]

natural = Word(nums)
parse_float = Regex(r'(\-)?\d+(\.)(\d*)?([eE][\-\+]\d+)?')
header = natural + LineEnd() + restOfLine
label = Word(alphas, max=2)
xyz = parse_float * 3

parse_atom = label.setResultsName("label") + xyz.setResultsName("xyz")
parser_xyz = Suppress(header) + OneOrMore(Group(parse_atom))
