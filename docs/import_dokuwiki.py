#!/usr/bin/env python3

import re
import sys
import subprocess


# ======================================================================================
def linker_old(match):
    txt, ref = match.groups()
    txt = txt.replace(r"\_", r"_")

    if ref[0].isdigit():
        link = f"https://dx.doi.org/{ref}"
    elif ref.split("/", 1)[0] == "CP2K_INPUT":
        ref = ref.replace("#list_", "/")
        ref = ref.replace("/", ".")
        link = f"#{ref}"
    elif ref.split("/", 1)[0] in {"FORCE_EVAL", "GLOBAL", "MOTION"}:
        ref = ref.replace("#list_", "/")
        ref = ref.replace("#", "/")
        ref = ref.replace("/", ".")
        link = f"#CP2K_INPUT.{ref}"
    elif ref.split("/", 1)[0] in {"tools", "src", "data"}:
        link = f"https://github.com/cp2k/cp2k/tree/master/{ref}"
    else:
        raise Exception(f"Unknown reference: {ref}")

    return f"[{txt}]({link})"


# ======================================================================================
def main_old():
    input_fn = sys.argv[1]
    cmd = ["pandoc", "--from=dokuwiki", "--to=markdown", input_fn]
    p = subprocess.run(cmd, capture_output=True, check=True)
    markdown = p.stdout.decode("utf8")
    markdown = re.sub(
        r"\[(.*?)\]\(https://www.google.com/search\?q=(.*?)&btnI=lucky\)",
        linker,
        markdown,
    )
    print(markdown)


# ======================================================================================
def convert_link(dokuwiki_link):
    inner = re.match("\[\[(.+)\]\]", dokuwiki_link.group(0)).group(1)
    parts = inner.split("|")
    target = parts[0].strip()
    label = parts[1].strip() if len(parts) > 1 else ""
    if target.startswith("doi>"):
        return f"[{label}](https://dx.doi.org/{target[4:]})"
    elif target.startswith("inp>"):
        path = target[4:].replace("/", ".").replace("#", ".")
        return f"[{label}](#CP2K_INPUT.{path})"
    elif target.startswith("https://"):
        return f"[{label}]({target})"
    else:
        print(f"Unknown target {target}")
        return f"[{label}](TODO:{target})"


# ======================================================================================
def convert_image(dokuwiki_image):
    inner = re.match("{{(.+)}}", dokuwiki_image.group(0)).group(1)
    parts = inner.split("|")
    target = parts[0].strip()
    caption = parts[1].strip() if len(parts) > 1 else ""

    filename = target.split(":", 1)[1].rsplit("?", 1)[0]
    width = target.rsplit("?", 1)[1] if "?" in target else ""

    if target.endswith(".zip"):
        return f"[{caption}](https://www.cp2k.org/_media/{target})"

    attributes = "align=center"
    if width:
        attributes += f" width={width}px"
    img_directive = f"![{caption}](images/{filename}){{{attributes}}}"

    return f"\n{img_directive}\n\n"

    # if not caption:
    #    return f"\n{img_directive}\n\n"
    #
    ## https://myst-parser.readthedocs.io/en/latest/syntax/optional.html#syntax-md-figures
    # xref = filename.rsplit(".", 1)[0] + "-figure"
    # output = [
    #    "",
    #    "",
    #    f":::{{figure-md}} {xref}",
    #    img_directive,
    #    "",
    #    caption,
    #    ":::",
    #    "",
    #    "",
    # ]
    # return "\n".join(output)


# ======================================================================================
def main():
    input_fn = sys.argv[1]
    assert input_fn.endswith(".dw")
    content = open(input_fn).read()
    content = re.sub(r"^====== ([^=]+) =+\s*$", r"# \1", content, flags=re.MULTILINE)
    content = re.sub(r"^===== ([^=]+) =+\s*$", r"## \1", content, flags=re.MULTILINE)
    content = re.sub(r"^==== ([^=]+) =+\s*$", r"### \1", content, flags=re.MULTILINE)
    content = re.sub(r"^=== ([^=]+) =+\s*$", r"#### \1", content, flags=re.MULTILINE)
    content = re.sub(r"^== ([^=]+) =+\s*$", r"##### \1", content, flags=re.MULTILINE)

    content = re.sub(r"^<code.*$", r"\n```", content, flags=re.MULTILINE)
    content = re.sub(r"^</code>\s*$", r"\n```", content, flags=re.MULTILINE)
    content = re.sub(r"^<file.*$", r"```\n", content, flags=re.MULTILINE)
    content = re.sub(r"^</file>\s*$", r"```\n", content, flags=re.MULTILINE)

    content = re.sub(r"^\\begin\{equation\}\s*$", r"$$", content, flags=re.MULTILINE)
    content = re.sub(r"^\\end\{equation\}\s*$", r"$$", content, flags=re.MULTILINE)

    content = re.sub(r"\[\[.+?\]\]", convert_link, content)
    content = re.sub(r"{{.+?}}", convert_image, content)

    output_fn = input_fn.replace(".dw", ".md")
    with open(output_fn, "w") as f:
        f.write(content)
    print(f"Wrote {output_fn}")


# ======================================================================================
main()

# EOF
