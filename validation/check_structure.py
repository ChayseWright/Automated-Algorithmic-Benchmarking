"""Small lexical/bracket/block audit. This is NOT the MATLAB parser/checkcode."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def clean(src):
    out = []
    i = 0
    while i < len(src):
        c = src[i]
        if c == "%":
            stop = src.find("\n", i)
            i = len(src) if stop < 0 else stop
            continue
        if c == "'":
            # All transpose operators in this package directly follow an operand.
            previous = src[i-1] if i else ""
            if previous and (previous.isalnum() or previous in ")]}."):
                i += 1
                continue
            i += 1
            while i < len(src):
                if src[i] == "'":
                    if i+1 < len(src) and src[i+1] == "'":
                        i += 2
                        continue
                    i += 1
                    break
                i += 1
            else:
                raise AssertionError("Unclosed single-quoted string")
            out.append(" STRING ")
            continue
        out.append(c)
        i += 1
    return "".join(out)


count = 0
for file in ROOT.rglob("*.m"):
    try:
        src = clean(file.read_text())
    except AssertionError as exc:
        raise AssertionError(str(file) + ': ' + str(exc)) from exc
    brackets, blocks = [], []
    for token in re.findall(r"[A-Za-z_]\w*|[()\[\]{}]", src):
        if token in ("(","[","{"):
            brackets.append(token)
        elif token in (")","]","}"):
            assert brackets and brackets.pop() == {")":"(","]":"[","}":"{"}[token], (file,token)
        elif not brackets:
            if token in ("function","if","for","while","switch","try","parfor","classdef"):
                blocks.append(token)
            elif token == "end":
                assert blocks, (file,"unmatched end")
                blocks.pop()
    assert not brackets and not blocks, (file,brackets,blocks)
    first = re.search(r"^function\s+(?:(?:\[[^\]]+\]|\w+)\s*=\s*)?(\w+)\s*\(",src)
    if first:
        assert first.group(1) == file.stem, (file,first.group(1))
    count += 1
print(f"Lexical/bracket/block audit passed for {count} MATLAB files (not a native syntax check).")
