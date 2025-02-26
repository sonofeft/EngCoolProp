import os, sys
import subprocess
import glob

try:
    vers_file = glob.glob( os.path.join( '*', '_version.py' ) )[0]
    print( vers_file )
except:
    print( "FAILED to find _version.py file." )
    sys.exit()

proc = subprocess.Popen(["git", "log", "--", vers_file], 
       stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()
s_out = out.decode('ascii')
#print("s_out output:", s_out)

lineL = []
sL = s_out.split( "\n" )

for line in sL:
    if line.startswith("commit"):
        lineL.append( line.strip() ) 

def get_v_from_s( s ):
    try:
        version = s.split()[2]
    except:
        version = ''
    return version.replace( "'","" ).replace('"','')

version_tagD = {} # key=version id: value=version number
for line in lineL:
    id = line.split()[-1]
    print( id )
    
    proc = subprocess.Popen(["git", "diff", id, "HEAD", "--", vers_file], 
           stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    diff_out = out.decode('ascii')
    #print("diff_out output:", diff_out)
    
    diffL = diff_out.split( "\n" )
    for s in diffL:
        if s.startswith( "+__version__ =" ):
            version_tagD["HEAD"] = get_v_from_s( s )
        elif s.startswith( "-__version__ =" ):
            version_tagD[ id ] = get_v_from_s( s )
            
print()
for id,v in version_tagD.items():
    print(id, v )
    
    # git tag -a v1.4 -m "my version 1.4" 9fceb02
    tag_cmd = 'git tag -a "%s" -m "version %s" %s' % (v, v, id)
    
    print( tag_cmd )
    # Un-Comment os.system in order to create tag commands
    # os.system( tag_cmd )
    
    
print( "="*66 )
print( "In order to push annotated tags to GitHub" )
#print( "git push origin --tags" )
print( "git push --follow-tags" )
print( "(NOTE: will require browser login)" )
