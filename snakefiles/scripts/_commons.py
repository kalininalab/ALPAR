from pathlib import Path
from typing import Any, Mapping, Generator, NoReturn, overload


@overload
def force_new_file(file: Path | str) -> Path: ...
@overload
def force_new_file(file: Any) -> NoReturn: ...

def force_new_file(file):
    """Enforce a filepath to be deleted.

    Parameters
    ----------
    file : Path | str | Any
        Path to file
    
    Returns
    -------
    Path
        Path to non-existent file.
    
    Raises
    ------
    ValueError
        If input value cannot be parsed as a Path object.
    """
    try:
        f = Path(file)
    except TypeError as e:
        raise ValueError(f"Could not parse path {file} to a Path object.") from e
    f.unlink(missing_ok=True)
    return f


def reverse_mapping[K, V](mapping: Mapping[K, V]) -> Mapping[V, K]:
    """Reverse a dictionary mapping.
    
    Parameters
    ----------
    mapping : Mapping[K, V]
        The mapping to reverse.
    
    Returns
    -------
    Mapping[V, K]
        The reversed mapping.
    
    Raises
    ------
    AssertionError
        If the mapping values are not unique.
    
    Examples
    --------
    >>> reverse_mapping({'a': 1, 'b': 2})
    {1: 'a', 2: 'b'}
    """
    assert len(mapping) == len(set(mapping.values())), "Mapping values are not unique."
    return {v: k for k, v in mapping.items()}


def zip_header_and_concat_content(file: Path, sep: str = '') -> Generator[tuple[str, str], None, None]:
    """Zip headers and concatenated content from a FASTA file.
    
    Parameters
    ----------
    file : Path
        The Path object to the FASTA file.
    sep : str, optional
        The separator to use when concatenating content lines, by default ''.
    
    Yields
    ------
    Generator[tuple[str, str], None, None]
        A generator yielding tuples of (header, concatenated content).
    """
    header = ''
    content = []
    with file.open('r', encoding='utf-8') as f:
        for line in iter(f):
            line = line.strip()
            if line.startswith('>'):
                if header:
                    yield header, sep.join(content)
                header = line
                content = []
            else:
                content.append(line)
        else:
            yield header, sep.join(content)


@overload
def write_to_file(file: tuple[Path | str, str]) -> None: ...
@overload
def write_to_file(file: Path | str, content: str = '') -> None: ...

def write_to_file(file, content = '') -> None:
    """Write content to a file.

    Makes sure the parent directory exists before writing.
    
    Parameters
    ----------
    file : Path | str
        The path of the file to write.
    content : str
        The content to write to the file.    
    """
    if isinstance(file, tuple):
        file, content = file
    f = Path(file)
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(content)
