from tests.test_imports import import_all_modules, import_lite_modules
from sr_amr.utils import is_tool_installed
from sr_amr import main

def test_imports():
    assert import_all_modules()
    
def test_imports_lite():
    assert import_lite_modules()

def test_tools():
    for tool in ["snippy", "prokka", "panaroo", "PanACoTA", "pyseer", "mashtree", "iqtree", "datasail"]:
        assert is_tool_installed(tool)

def test_tools_lite():
    for tool in ["snippy"]:
        assert is_tool_installed(tool)

def test_strain_file_creation():
    

