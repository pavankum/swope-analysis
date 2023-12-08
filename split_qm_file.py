with open('qm-all-annotated-full.sdf','r') as readfile:
    lines = readfile.readlines()[:400]

    newfile = []
    for i,line in enumerate(lines):
        # print(line)
        if line.find('$$$$') != -1:
            newfile.append('> <final_energy>\n') # need energy listed as "final_energy"
            newfile.append(energy+'\n')
            newfile.append(line)

            # write previous file
            # print(title)
            with open('b3lyp-d3bj_dzvp/'+title+'.sdf','w') as writefile:
                writefile.write(title+'\n') # Needs title to be the compound ID
                for line2 in newfile[1:]:
                    writefile.write(line2)
            newfile = []

        else:
            newfile.append(line)
            if line.find('Compound Name') != -1:
                nameline = lines[i+1].strip('\n')
                splitline =nameline.split('-')
                conf_num = splitline[-1].split('.')[0]

                name_noconf = '-'.join(splitline[:-1])
                if len(conf_num) <2:
                    conf_num = '0'+conf_num
                title = name_noconf + '-' + conf_num# + '.sdf'
            if line.find('Energy QCArchive') != -1:
                energy = lines[i+1]
