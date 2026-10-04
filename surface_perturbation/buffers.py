"""Share numerical buffers for exact duplicate files, retaining every run identity."""


def repeated_field(pair,buffer,compute):
    key=(pair.get('wt_sha256'),pair.get('mut_sha256'),pair.get('chain'))
    if not all(key[:2]):return compute(),False
    if key in buffer:return buffer[key],True
    result=compute()
    buffer[key]=result
    return result,False
