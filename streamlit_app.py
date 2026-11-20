  import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(page_title="Spill-Proof QR", layout="wide")

# Read your teammate's HTML and Python files directly from the repo root
with open("index.html", "r", encoding="utf-8") as f:
    html_code = f.read()

with open("qr_engine.py", "r", encoding="utf-8") as f:
    py_code = f.read()

# Inject qr_engine.py inline so Pyodide doesn't fail on relative file fetches inside Streamlit's iframe
escaped_py = py_code.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")
html_code = html_code.replace(
    "let srcUrl = 'qr_engine.py', res = await fetch(srcUrl).catch(() => null);\n        if (!res || !res.ok) { srcUrl = 'qr_engine.txt'; res = await fetch(srcUrl); }\n        const src = await res.text();",
    f"const src = `{escaped_py}`;"
)

components.html(html_code, height=1400, scrolling=True)
