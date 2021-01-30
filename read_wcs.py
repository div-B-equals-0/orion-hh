import toml

def read_wcs(wcsfile: str) -> dict: 
    "Read WCS (key,value) pairs from DS9's *.wcs files"
    with open(wcsfile) as stream:
        return toml.loads(_fixups(stream.read()))


def _fixups(s: str) -> str:
    s = s.replace("0.\n", "0.0\n")
    return s


if __name__ == "__main__":
    data = read_wcs("mosaic-1996-HH204-align-robberto.wcs")
    print(data)
