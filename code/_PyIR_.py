from crowelab_pyir import PyIR

def main():
    FILE = 'tmp/sequence.fasta'
    ##save as tsv file
    #pyir = PyIR(query=FILE, args=['--outfmt', 'tsv','-o','tmp/output'])    #equal to shell: pyir tmp/sequence.fasta --outfmt tsv -o tmp/output
    ##pasta as dictionary
    pyir = PyIR(query=FILE, args=['--outfmt', 'dict'])
    result = pyir.run()
    print(result.items())
    print(len(result))
if __name__ == "__main__":
    main()