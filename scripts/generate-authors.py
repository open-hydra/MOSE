#!/usr/bin/env python3

from pathlib import Path
import subprocess
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
CONFIG = SCRIPT_DIR / "contributors.yml"
OUTPUT = REPO_ROOT / "AUTHORS.md"

with open(CONFIG, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f)

def normalize_people(items):
    people = []

    for item in items:
        if isinstance(item, str):
            people.append({"name": item, "role": ""})
            continue

        if isinstance(item, dict) and "name" in item:
            people.append(
                {
                    "name": item["name"],
                    "role": item.get("role", ""),
                }
            )

    return people

original_authors = normalize_people(cfg.get("original_authors", []))
maintainers = normalize_people(
    cfg.get("present_maintainers", cfg.get("maintainers", []))
)
alumni = cfg.get("alumni", [])
exclude = set(cfg.get("exclude", []))

result = subprocess.check_output(
    ["git", "shortlog", "-sne", "--all"],
    text=True,
)

contributors = []

for line in result.splitlines():
    line = line.strip()

    if not line:
        continue

    # Match previous pipeline filtering that excluded bot-generated identities.
    if "dependabot" in line.lower() or "github-actions" in line.lower():
        continue

    _, identity = line.split("\t", 1)

    name = identity.split("<")[0].strip()

    if name in exclude:
        continue

    contributors.append(name)

maintainer_names = {m["name"] for m in maintainers}
original_author_names = {a["name"] for a in original_authors}
alumni_names = set(alumni)

contributors = sorted(
    set(contributors)
    - maintainer_names
    - original_author_names
    - alumni_names
)

with open(OUTPUT, "w", encoding="utf-8") as f:
    f.write("# Authors\n\n")

    if original_authors:
        f.write("## Original Authors\n\n")

        for a in original_authors:
            if a["role"]:
                f.write(f'- {a["name"]} - {a["role"]}\n')
            else:
                f.write(f'- {a["name"]}\n')

        f.write("\n")

    f.write("## Present Maintainers\n\n")

    for m in maintainers:
        if m["role"]:
            f.write(f'- {m["name"]} - {m["role"]}\n')
        else:
            f.write(f'- {m["name"]}\n')

    f.write("\n## Contributors\n\n")
    f.write("Automatically sourced from git history.\n\n")

    for c in contributors:
        f.write(f"- {c}\n")

    if alumni:
        f.write("\n## Alumni\n\n")

        for a in alumni:
            f.write(f"- {a}\n")

print(f"Generated {OUTPUT}")