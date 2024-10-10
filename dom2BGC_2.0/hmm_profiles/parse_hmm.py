import argparse

def commandLineParser():
    '''
    Parses input and holds default values.
    Returns a dictionary with all the parameters for the wrapper
    '''
    parser = argparse.ArgumentParser(description='Parses hmmsearch output of natural products associated domains.')
    #General argnuments to specify where to put/find files
    parser.add_argument('-i', '--hmmSearchOutput', type=str, required=True, help='full path to the hmmsearch output file')
    parser.add_argument('-o', '--output', type=str, required=True, help='full path to store the parsed output file')
    parser.add_argument('-pos', '--startEndPosition', required=False, default='169,237', help='comma separated list of start and end position for the hmm profile, default="169,237"')
    return vars(parser.parse_args())

def parse_alignment(alignment, outfile, S, E):
    alignment = alignment.split('Domain annotation for each sequence (and alignments):\n')[1]
    alignment = alignment.split('\n\n\n')[0].split('>> ')[1:]
    amp_domains = ''
    for block in alignment:
        number_hits = block.split('\n\n')[0]
        lines = number_hits.split('\n')
        domains = block.split('== domain')[1:]
        for x in range(3, len(lines)):
            line = [y for y in lines[x].split(' ') if y]
            start, finish = int(line[6]), int(line[7])
            print (start, finish)
            if start <= S+5 and finish >= E-5:
                seq = get_alignment(domains[x-3])
                seq = seq[S:E]
                print (seq)
                if len(seq) >= E-S:
                    amp_domains += '>{}_{}\n{}\n'.format(block.split(' ')[0], x-2, seq)
    out_txt = open(outfile, 'w')
    out_txt.write(amp_domains)
    out_txt.close()

def get_alignment(domain):
    lines = domain.split('\n')
    l = [x for x in lines[2].split(' ') if x]
    dom_seq = '-'*(int(l[1])-1)
    for y in range(4, len(lines), 6):
        dom_seq+=[x for x in lines[y].split(' ') if x][2]
    dom_seq = ''.join([x for x in dom_seq if not x.islower()])
    return dom_seq

if __name__ == '__main__':
    argOptions = commandLineParser()
    infile = open(argOptions['hmmSearchOutput']).read()
    outfile = argOptions['output']
    S = int(argOptions['startEndPosition'].split(',')[0])
    E = int(argOptions['startEndPosition'].split(',')[1])
    parse_alignment(infile, outfile, S, E)