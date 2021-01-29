import os
import sys
import subprocess
from pathlib import Path

venv_path = Path(os.environ['VIRTUAL_ENV'])
notebook_path = venv_path.parent
jupyter_cmd = 'jupyter' # should use the venv's jupyter, right

def disable_use_redirect_file():
    jupyter_config = Path("~/.jupyter/jupyter_notebook_config.py").expanduser()
    if not jupyter_config.exists():
        subprocess.check_call([jupyter_cmd, "notebook", "--generate-config"])
    fix = "\\nc.NotebookApp.use_redirect_file = False\\n"
    if not fix in jupyter_config.read_text():
        with open(jupyter_config, "a+") as op:
            op.write_text(fix)
    print("Disabled file based jupyter redirect")

def place_shortcut_on_desktop():
    if sys.platform == 'Windows':
        cmd = 'PowerShell -NoProfile --Command "Write-Host([Environment]::GetFolderPath(\'Desktop\'))"'
        desktop_folder = Path(subprocess.check_output(cmd, shell=True).decode('utf-8'))
        bin = 'Scripts'
    else:
        desktop_folder = Path("~/Desktop").expanduser()
        bin = 'bin'
    target = desktop_folder / ("jupyter notebook " + notebook_path.name + '.py')
    jcmd = venv_path / bin / 'jupyter'
    target.write_text(f"""#!/usr/bin/env python3
import subprocess
subprocess.call(["{jcmd}", 'notebook'], cwd="{notebook_path}")
""")
    if sys.platform != 'Windows':
        os.chmod(target, 0o755)
    print("Placed jupyter shortcut on desktop")

def main():
    if sys.platform == 'Windows':
        disable_use_redirect_file()
    place_shortcut_on_desktop()


