""" Cython implementation of a CDB reader """
from libc.stdio cimport fopen, FILE, getline, fclose, sscanf, fscanf, fread, fseek, SEEK_CUR
import numpy as np
cimport numpy as np
from libc.stdlib cimport atoi, atof
from libc.string cimport strncpy, strcmp



cimport cython
@cython.boundscheck(False) # turn of bounds-checking for entire function
@cython.wraparound(False)
def Read(filename):
    filename_byte_string = filename.encode("UTF-8")
    cdef char* fname = filename_byte_string
    
    # Check file exists
    cdef FILE* cfile
    cfile = fopen(fname, "rb")
    if cfile == NULL:
        raise Exception("No such file or directory: '%s'" % filename)
 
    # Define variables
    cdef char * line = NULL
    cdef size_t l = 0
    cdef ssize_t read
    cdef int[5] blocksz
    cdef int i, j, nlines
    cdef int tempint
    cdef int nnodes, linelen, isz
    cdef float tempflt

    # Size temp char array
    cdef char tempstr[100]
    
    # Get element types
    elem_type = []
    rnum = []
    rdat = []
    while True:
        getline(&line, &l, cfile)
        
        # Record element types
        if 'E' == line[0]:
            if 'ET' in line:
                elem_type.append([int(line[3:line.find(',', 5)]),     # element number
                                  int(line[line.find(',', 5) + 1:])]) # element type        
        
        if 'R' == line[0]:
            if 'RLBLOCK' in line:
                # Get number of sets
                ist = line.find(',') + 1
                ien = line[ist:].find(',') + ist
                nset = int(line[ist:ien])
            
                # Skip Format1 and Format2 (always 2i8,6g16.9 and 7g16.9)
                getline(&line, &l, cfile)
                getline(&line, &l, cfile)
                
                # Read data
                c_set = 0
                while True:
                    getline(&line, &l, cfile)
                    rcon = [] # real constants
                    
                    c_set += 1
                    if c_set > nset:
                        break
                    
                    # Get real constant number
                    rnum.append(int(line[:8]))
                    
                    # Number of constants
                    ncon = int(line[8:16])
                    
                    # Get constant data
                    if ncon > 6: # if multiple lines
                        for i in range(6):
                            rcon.append(float(line[16 + 16*i:32 + 16*i]))
                            ncon -= 1
                            
                        # advance line
                        getline(&line, &l, cfile)
                         
                        # read next line
                        while True:
                            if ncon > 7:
                                for i in range(7):
                                    rcon.append(float(line[16*i:16*(i + 1)]))
                                    ncon -= 1
                                # advance
                                getline(&line, &l, cfile)
                                    
                            else:
                                for i in range(ncon):
                                    rcon.append(float(line[16*i:16 + 16*i]))  
                                    
                                break
                            
                    # If only one in constant data
                    else:
                        for i in range(ncon):
                            rcon.append(float(line[16 + 16*i:32 + 16*i]))   
            
                    rdat.append(rcon)
        
        if 'N' == line[0]: # Test is faster than next line
        
            # if line contains the start of the node block
            if 'NBLOCK' in line:
                # Get size of NBLOCK
                nnodes = int(line[line.rfind(',') + 1:])
                
                # Get format of NBLOCk
                getline(&line, &l, cfile)
                d_size, f_size, nfld = GetBlockFormat(line)
                break

            
    # Size node arrays
    cdef int [::1] nnum = np.empty(nnodes, dtype=np.int32)
    cdef double [:, ::1] nodes = np.empty((nnodes, 6))
    nodes[:] = 0.0
    
    # Create the interger string and place the null character
    cdef int intsize = d_size
    cdef char intstr[100]
    intstr[intsize + 1] = '\0'
    
    # Create floating point number string and place the null character
    cdef int fsize = f_size
    cdef char floatstr[100]
    floatstr[fsize + 1] = '\0'

    # Loop through all nodes and store number of fields
    cdef int nfields = nfld
    for i in range(nnodes):
        # Read node number from field 1
        fread(intstr, 1, intsize, cfile)
        nnum[i] = atoi(intstr)
        
        # skip fields 2 and 3
        fseek(cfile, intsize*2, SEEK_CUR)
        
        # Read fields 4-9 (or 4-6)
        for j in range(nfields):
            fread(floatstr, 1, fsize, cfile)
            if '\r' == floatstr[0]: # if dos EOL
                # Seek backwards as we've gone past the EOL
                fseek(cfile, -(fsize - 2), SEEK_CUR)
                break
            
            elif '\n' == floatstr[0]:
                # Seek backwards as we've gone past the EOL
                fseek(cfile, -(fsize - 1), SEEK_CUR)
                break
            
            else:
                # Otherwise, convert and store
                nodes[i, j] = atof(floatstr)
            
    
        # If all fields have been read
        if j == nfields - 1:
            fread(tempstr, 1, 1, cfile)
            if '\r' == tempstr[0]: # if dos EOL
                # read the next eol character and continue
                fseek(cfile, 1, SEEK_CUR)


    ############### EBLOCK ###############
    # Seek to the start of the element data
    while True:

        getline(&line, &l, cfile)
#        if read == -1:
#            break
        
        if 'E' == line[0]:
        
            # if line contains the start of the node block
            if 'EBLOCK' in line:
                # Get size of EBLOCK
                nlines = int(line[line.rfind(',') + 1:])
                
                # Get interger block size
                getline(&line, &l, cfile)
                isz = int(line[line.find('i') + 1:line.find(')')])
                break
            

            
    # Initialize element data array.  Use number of lines as nelem is unknown
    cdef int [:, ::1] elem = np.empty((nlines, 20), dtype=np.int32, order='C')
    cdef int [::1] etype = np.empty(nlines, dtype=np.int32, order='C')
    cdef int [::1] elemnum = np.empty(nlines, dtype=np.int32, order='C')
    cdef int [::1] e_rcon = np.empty(nlines, dtype=np.int32, order='C')
    cdef int nnode, nelem

    # Null element is -1
    elem[:] = -1

    i = 0 # init counter
    while True:
    
        # Check if at end of EBLOCK
        fscanf(cfile, '%d', &tempint)
        if tempint == -1:
            break
        
        # Field 2: Read element type
        fscanf(cfile, '%d', &tempint)
        etype[i] = tempint
    
        # Field 3: Read real constant
        fscanf(cfile, '%d', &tempint)
        e_rcon[i] = tempint

        # Skip Fields 4 - 8 and store 9, the number of nodes    
        for c in range(6):
            fscanf(cfile, '%d', &nnode)

        # Store element number
        for c in range(2):
            fscanf(cfile, '%d', &tempint)
        elemnum[i] = tempint        
            
        # Read nodes in element
        for c in range(nnode):
            fscanf(cfile, '%d', &tempint)
            elem[i, c] = tempint
            
        # next element
        i += 1
        
    # Store actual number of elements
    nelem = i


    # Get node components
    cdef int ncomp
    cdef int [::1] component
    cdef int nblock
        
    # Store node compondents
#    node_comps = {'FIXED_NODES': [], 'NO_CHECK': []} # always include these
    node_comps = {}
    while True:        
        read = getline(&line, &l, cfile)
        if read == -1:
            break
    
        if 'C' == line[0]:
            if 'CMBLOCK' and 'NODE' in line:

                # Get Component name
                ind1 = line.find(',') + 1
                ind2 = line.find(',', ind1)
                comname = line[ind1:ind2]

                # Get number of items
                ncomp = int(line[line.rfind(',') + 1:line.find('!')])
                component = np.empty(ncomp, dtype=np.int32, order='C')
                
                # Get interger size
                getline(&line, &l, cfile)
                isz = int(line[line.find('i') + 1:line.find(')')])
                tempstr[isz] = '\0'
                
                # Number of intergers per line
                nblock = int(line[line.find('(') + 1:line.find('i')])
                
                # Extract nodes
                for i in xrange(ncomp):
                    
                    # Read new line if at the end of the line
                    if i%nblock == 0:
                        getline(&line, &l, cfile)
                    
                    strncpy(tempstr, line + isz*(i%nblock), isz)
                    component[i] = atoi(tempstr)

                # Convert component to array and store
                node_comps[comname] = ComponentInterperter(component)
                      
    # Close file
    fclose(cfile)

    return {'rnum': np.asarray(rnum),
            'rdat': np.asarray(rdat),
            'ekey': np.asarray(elem_type),
            'nnum': np.asarray(nnum),
            'nodes': np.asarray(nodes),
            'enum': np.asarray(elemnum[:nelem]),
            'elem': np.array(elem[:nelem]),
            'etype': np.asarray(etype[:nelem]),
            'e_rcon': np.asarray(e_rcon),
            'node_comps': node_comps}
     
    
def GetBlockFormat(string):
    """ Get node block format """
    
    # Digit Size
    d_size = int(string[string.find('i') + 1:string.find(',')])
    f_size = int(string[string.find('e') + 1:string.find('.')])
    nfields = int(string[string.find(',') + 1:string.find('e')])

    return d_size, f_size, nfields

            
def ComponentInterperter(component):
    """
    If a node is negative, it is describing a list from the previous node.  This is ANSYS's way of 
    saving file size when writing components.
    
    """
    
    f_new = []
    for i in range(len(component)):
        if component[i] > 0: # Append if positive
            f_new.append(component[i])
        else: # otherwise, append list
            f_new.append(range(abs(component[i - 1]) + 1, abs(component[i]) + 1))
    
    return np.hstack(f_new).astype(np.int32)
