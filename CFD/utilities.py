from decimal import Decimal

def importobj(file,dl):
    with open(file) as f:
        l = f.readlines()
    p=[]
    for _ in l:
        if _[0]=='v' and _[1]==' ':
            p.append(_.split()[1:])
    for i,v in enumerate(p):
        p[i]=(int(Decimal(v[0])/Decimal(dl)),int(Decimal(v[1])/Decimal(dl)))
    return p

def status(it,re,t,cli,mode):
    if mode=='run':
        if cli:
            print(f'iteration:       {it}')
            print(f'U residual:      {re:.5f}')
            print(f'elapsed time:    {t//60:.0f}m {t%60:.0f}s   ')
            print('\033[1A\033[1A\033[1A',end='\x1b[2K')
        else:
            print((f'iteration: {it} | U residual: {re:.5f} | elapsed time: {t//60:.0f}m {t%60:.0f}s   '),end='\r')     
    elif mode=='start':
        if cli: 
            print('\033[?25l',end='')
    elif mode=='end':
        if cli: 
            print('\033[?25h',end='')
            print(f'iteration:       {it}')
            print(f'U residual:      {re:.5f}')
            print(f'elapsed time:    {t//60:.0f}m {t%60:.0f}s   ')
        else:
            print((f'iteration: {it} | U residual: {re:.5f} | elapsed time: {t//60:.0f}m {t%60:.0f}s   '))

def in_cli():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell': return False # Jupyter
        elif shell == 'TerminalInteractiveShell': return True
        else: return True
    except NameError:
        return True