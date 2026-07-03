"""Render README.md, learning-plan.md, and each module's concepts.md to
styled .html twins with MathJax support.

GitHub Pages serves raw .md files as plain text (content-type:
text/markdown), so this script produces a rendered .html twin next to each
source .md for the site to link to instead. The .md files remain the
source of truth. Run this script after any edit to README.md,
learning-plan.md, any module's concepts.md, or this template:

    python3 tools/render_md_pages.py
"""

import markdown
import re
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXT_REPO_PAGES = "https://fhoces.github.io/experimentation-refresher/"

TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    max-width: 760px;
    margin: 3em auto;
    padding: 0 1.2em;
    color: #222;
    line-height: 1.6;
  }}
  h1 {{ border-bottom: 1px solid #ddd; padding-bottom: 0.3em; }}
  h2 {{ margin-top: 1.8em; }}
  h3 {{ margin-top: 1.4em; }}
  a {{ color: #0366d6; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1.2em 0; }}
  th, td {{ padding: 0.5em 0.75em; border: 1px solid #e1e4e8; text-align: left; }}
  th {{ background: #f6f8fa; }}
  code {{ background: #f6f8fa; padding: 0.15em 0.4em; border-radius: 3px; font-size: 0.9em; }}
  pre {{ background: #f6f8fa; padding: 1em; border-radius: 5px; overflow-x: auto; }}
  pre code {{ background: none; padding: 0; }}
  blockquote {{ border-left: 4px solid #ddd; margin-left: 0; padding-left: 1em; color: #555; }}
  .site-nav {{ font-size: 0.85em; margin-bottom: 1.5em; color: #888; }}
</style>
<script>
MathJax = {{
  tex: {{ inlineMath: [['$', '$']], displayMath: [['$$', '$$']] }}
}};
</script>
<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
</head>
<body>
{nav}
{body}
</body>
</html>
"""

ROOT_NAV = '<nav class="site-nav"><a href="index.html">Course home</a></nav>'
MODULE_NAV = (
    '<nav class="site-nav"><a href="../index.html">Course home</a> '
    '&middot; <a href="slides.html">Slides for this module</a></nav>'
)


def rewrite_links(src_text):
    # GitHub repo link for the sibling course -> the sibling's own Pages site
    src_text = src_text.replace(
        "(https://github.com/fhoces/experimentation-refresher)", f"({EXT_REPO_PAGES})"
    )
    # ../experimentation-refresher/ -> absolute Pages URL of the sibling repo
    src_text = src_text.replace("(../experimentation-refresher/)", f"({EXT_REPO_PAGES})")
    # learning-plan.md -> learning-plan.html (only at repo root context)
    src_text = re.sub(r'\]\(learning-plan\.md\)', '](learning-plan.html)', src_text)
    # module-XX/concepts.md -> module-XX/concepts.html
    src_text = re.sub(r'\]\((module-\d+)/concepts\.md\)', r'](\1/concepts.html)', src_text)
    return src_text


def render(md_path, out_path, title, nav):
    with open(md_path) as f:
        text = f.read()
    text = rewrite_links(text)
    body = markdown.markdown(text, extensions=["extra", "sane_lists", "toc"])
    html = TEMPLATE.format(title=title, nav=nav, body=body)
    with open(out_path, "w") as f:
        f.write(html)
    print(f"wrote {out_path} ({len(html)} bytes)")


# README.md
render(f"{REPO}/README.md", f"{REPO}/README.html", "Causal Inference Beyond A/B Tests", ROOT_NAV)

# learning-plan.md
render(
    f"{REPO}/learning-plan.md",
    f"{REPO}/learning-plan.html",
    "Learning Plan: Causal Inference Beyond A/B Tests",
    ROOT_NAV,
)

# each module's concepts.md
titles = {
    "module-01": "Module 1: TWFE Diagnosed",
    "module-02": "Module 2: Heterogeneity-Robust DiD",
    "module-03": "Module 3: Honest DiD",
    "module-04": "Module 4: Synthetic Control",
    "module-05": "Module 5: Synthetic DiD",
    "module-06": "Module 6: Causal Forest",
    "module-07": "Module 7: Policy Learning",
    "module-08": "Module 8: Matrix Completion",
}
for mod, title in titles.items():
    render(f"{REPO}/{mod}/concepts.md", f"{REPO}/{mod}/concepts.html", f"{title} - Concepts", MODULE_NAV)
