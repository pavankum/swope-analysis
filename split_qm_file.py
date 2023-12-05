with open('qm-all-annotated.sdf','r') as readfile:
    lines = readfile.readlines()

    newfile = []
    for i,line in enumerate(lines):
        # print(line)
        if line.find('$$$$') != -1:
            # print(newfile)
            # write previous file
            print(title)
            with open('QM_files/'+title+'.sdf','w') as writefile:
                for line2 in newfile:
                    writefile.write(line2)
            newfile = []

        else:
            newfile.append(line)
            if line.find('Compound Name') != -1:
                title = lines[i+1].strip('\n')

