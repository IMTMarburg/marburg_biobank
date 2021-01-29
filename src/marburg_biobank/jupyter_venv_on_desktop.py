import os
import sys
import subprocess
from pathlib import Path

venv_path = Path(sys.prefix)
notebook_path = venv_path.parent
is_windows = sys.platform.lower() in ("win32", "windows")
if is_windows:
    jupyter_cmd = (
        venv_path / "Scripts" / "jupyter.exe"
    )  # should use the venv's jupyter, right
else:
    jupyter_cmd = venv_path / "bin" / "jupyter"  # should use the venv's jupyter, right


def disable_use_redirect_file():
    jupyter_config = Path("~/.jupyter/jupyter_notebook_config.py").expanduser()
    if not jupyter_config.exists():
        subprocess.check_call([jupyter_cmd, "notebook", "--generate-config"])
    fix = "\nc.NotebookApp.use_redirect_file = False\n"
    if not fix in jupyter_config.read_text():
        with open(jupyter_config, "a+") as op:
            op.write(fix)
    print("Disabled file based jupyter redirect")


def place_shortcut_on_desktop():
    if is_windows:
        cmd = "PowerShell -NoProfile -Command \"Write-Host([Environment]::GetFolderPath('Desktop'))\""
        desktop_folder = Path(subprocess.check_output(cmd, shell=True).decode("utf-8").strip())
    else:
        desktop_folder = Path("~/Desktop").expanduser()
    target = desktop_folder / ("jupyter notebook " + notebook_path.name + ".py")
    target.write_text(
        f"""#!/usr/bin/env python3
import subprocess
subprocess.call([r"{jupyter_cmd}", 'notebook'], cwd=r"{notebook_path}")
"""
    )
    if is_windows:
        os.chmod(target, 0o755)
    print("Placed jupyter shortcut on desktop")


def main():
    if is_windows:
        disable_use_redirect_file()
    place_shortcut_on_desktop()
