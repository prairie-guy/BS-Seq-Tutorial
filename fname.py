# C. Bryan Daniels
# 1/1/24


def fname(*args):
    """
    usage: fname(path,ext,stem) | fname(path,ext)
    Returns a str with {sample} wildcard interpolated into it
    """
    if len(args) == 2:
        path,ext = args
        return f"{path}/{{sample}}.{ext}"
    elif len(args) == 3:
        path,ext,stem = args
        return f"{path}/{{sample}}_{stem}.{ext}"
    else:
        print("usage: fname(path,ext,stem) | fname(path,ext)")

def expand_fnames(*args):
    """
    usage: expand_fnames(path,ext,stem,samples) | expand_fnames(path,ext,samples)
    Returns a list of strings with {sample} wildcard interpolated into it
    """
    if len(args) == 3:
        path,ext,samples = args
        return [f"{path}/{sample}.{ext}" for sample in samples]
    elif len(args) == 4:
        path,ext,stem,samples = args
        return [f"{path}/{sample}_{stem}.{ext}" for sample in samples]
    else:
        print("usage: expand_fnames(path,ext,stem,samples) | expand_fnames(path,ext,samples)")
