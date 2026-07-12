import os

from sr_amr.utils import conda_env_wrapper


def test_conda_env_wrapper_accepts_multiple_args_and_kwargs(monkeypatch):
    monkeypatch.setenv("ALPAR_SUBPROCESS", "1")

    def inner(*args, **kwargs):
        return args, kwargs

    wrapped = conda_env_wrapper("alpar-ml")(inner)
    result = wrapped("model", "table", output_dir="out")

    assert result == (("model", "table"), {"output_dir": "out"})
