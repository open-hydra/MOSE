#!/usr/bin/env python3
"""Set or add keys of a MOSE input.ini: set_ini.py FILE 'Section|key|value' ... (a new key goes right after the last
key = value line of its section). The file always ends with a newline."""
import re
import sys


def main():
    path = sys.argv[1]
    lines = open(path, newline='').read().split('\n')
    keyre = re.compile(r'^\s*[A-Za-z0-9_.+-]+\s*=')
    for spec in sys.argv[2:]:
        sec, key, val = spec.split('|')
        hdr = [i for i, ln in enumerate(lines) if ln.strip().startswith('[')]
        h = next((i for i in hdr if lines[i].strip() == '[' + sec + ']'), None)
        if h is None:
            lines += ['', '[' + sec + ']', key + ' = ' + val]
            continue
        end = next((i for i in hdr if i > h), len(lines))
        rep = [i for i in range(h + 1, end) if re.match(r'^\s*' + re.escape(key) + r'\s*=', lines[i])]
        if rep:
            lines[rep[0]] = key + ' = ' + val
        else:
            last = max([i for i in range(h + 1, end) if keyre.match(lines[i])] or [h])
            lines.insert(last + 1, key + ' = ' + val)
    text = '\n'.join(lines)
    if not text.endswith('\n'):
        text += '\n'
    open(path, 'w', newline='').write(text)


if __name__ == '__main__':
    main()
