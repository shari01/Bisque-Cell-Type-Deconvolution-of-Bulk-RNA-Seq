import subprocess
def test_cli_help():
    out = subprocess.run(["bisque-deconv", "-h"], capture_output=True, text=True)
    assert out.returncode == 0
    assert "usage" in out.stdout.lower()
