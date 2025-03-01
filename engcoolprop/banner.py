
        
def floatDammit( val=0.0 ):
    val = str( val )
    v = val.strip()
    if v.lower() in ['nan','inf','-inf']:
        if v.lower() == '-inf':
            v = -1.0E99
        elif v.lower() == 'inf':
            v = 1.0E99
        else:
            v = 0.0
    try:
        v = float(v)
    except:
        v = 0.0
    return v

def is_float( fval ):
    if type(fval)==type(11.11):
        return 1
    if type(fval)==type('string'):
        if fval.find('.')==-1:
            return 0
    try:
        if float(fval) > -float("inf"):
            return 1
    except:
        #print 'Failed is_float',fval,type(fval)
        pass
    return 0

def get_banner_chars( outlineChar=None ):
    
    
    if outlineChar is None:# printable characters
        ul = '+'; ur = '+'; ll = '+'; lr = '+'; vert = '|'; hor = '-'
        vl = vr = vert
        ht = hb = hor
        um = lm = '+'
    else:
        # will be all outlineChar if not printable characters from above
        ul = ur = ll =  lr = vr = vl = hb = ht = vert = hor = um = lm = outlineChar
        
    return ul, ur, ll,  lr, vr, vl, hb, ht, vert, hor, um, lm


def banner(s, outlineChar=None, leftMargin=0, just='left'):
    
    ul, ur, ll,  lr, vr, vl, hb, ht, vert, hor, um, lm = get_banner_chars( outlineChar )
    
    sL = s.split('\n')
    L = 0
    for s in sL:
        ss = '  ' + s + '  '
        L = max(L, len(ss))
    
    top = leftMargin*' ' + ul + L*ht + ur
    bot = leftMargin*' ' + ll + L*hb + lr
    
    midL = []
    for s in sL:
        # if just=='center':
        #     midL.append( leftMargin*' ' + vl + '%s'%(s.center(L),) + vr )
        # elif just=='right':
        #     midL.append( leftMargin*' ' + vl + '%s'%(s.rjust(L),) + vr )
        # else:
        #     midL.append( leftMargin*' ' + vl + '%s'%(s.ljust(L),) + vr )

        midL.append( leftMargin*' ' + vl + '%s'%(s.center(L),) + vr )
        
    
    pad = ''
    if just=='center':
        Lpad = (80 - len(top)) // 2
        pad = ' '*Lpad
    elif just=='right':
        Lpad = 79 - len(top)
        pad = ' '*Lpad

    
    print(pad+top)
    for mid in midL:
        print(pad+mid)
    print(pad+bot)

def getStr( val ):
    if type(val) == type('str'):
        return val
    if is_float(val):
        val = floatDammit( val )
        #return ('%f' % val).rstrip('0').rstrip('.')
        return '%g'%val 
    else:
        return str( val )

def show_table( titleL=None, LOL=None, outlineChar=None, header_sep=True ):
    '''Number of columns = len(titleL)
       Lists of values in list of lists LOL
    '''
    #ul = chr(218); ur = chr(191); ll = chr(192); lr = chr(217); vert = chr(179); hor = chr(196)
    #um = chr(194); lm = chr(193)
    #vl = vr = vert
    #ht = hb = hor
    ul, ur, ll,  lr, vr, vl, hb, ht, vert, hor, um, lm = get_banner_chars( outlineChar )
    
    
    Ncol = len(titleL)
    lenL = [ len(str(t)) for t in titleL ] # allocation length for each column 
    # lenL.append( len(titleL[j]) ) # init to title Length
    for j in range( Ncol ):
        for L in LOL:
            val = L[j]
            lstr = len( getStr(val) )
            if lstr > lenL[j]: # set lenL[j] to max len in column
                lenL[j] = lstr
            

    #print 'lenL =',lenL
    # print top line 
    sL = [ul]
    for i in range(Ncol):
        sL.append( hor*(lenL[i]+2) )
        if i < Ncol-1:
            sL.append( um )
    sL.append(ur)
    print(''.join(sL))
    
    # print titles
    sL = []
    for i,title in enumerate(titleL):
        sL.append( title.center( lenL[i]+2 ) )
    print(vert + vert.join(sL) + vert)

    if header_sep:
        # print header separation line
        sL = [ul]
        for i in range(Ncol):
            sL.append( hor*(lenL[i]+2) )
            if i < Ncol-1:
                sL.append( um )
        sL.append(ur)
        print(''.join(sL))

    # print content
    for i in range( len(LOL) ):
        sL = []
        for j in range( Ncol ):
            sL.append( getStr(LOL[i][j]).center( lenL[j]+2 ) )
        print(vert + vert.join(sL) + vert)

    # print bottom line 
    sL = [ll]
    for i in range(Ncol):
        sL.append( hor*(lenL[i]+2) )
        if i < Ncol-1:
            sL.append( lm )
    sL.append(lr)
    print(''.join(sL))
    
def dev_tests():
        
    for i in range(170, 241, 10):
        for j in range(10):
            print(i+j,chr(i+j),' ', end=' ')
        print()
    banner('Default Banner')
    banner('outlineChar = *\njust = "left"', outlineChar='*')

    banner('outlineChar = None\njust = "center"', outlineChar=None, leftMargin=0, just='center')
    banner('outlineChar = #\njust = "right"', outlineChar='#', leftMargin=0, just='right')
    
    # --------------
    # for i in range(6,8):
    # use either default or i=0 or 7
    for c in [None,  '|', '.']:
        show_table( ['col1 outlineChar=%s'%c,'second Col','col 3'], [[1,2,3.333],['a','seven',123.4]], outlineChar=c )

    print()
    for v in ['nan','inf','-inf', 123, 99.0, 'hi']:
        print( 'v=%-4s'%v,'  getstr=%-4s'%getStr(v), 'floatDammit=',floatDammit(v))
    


if __name__=="__main__":
    dev_tests()    
    