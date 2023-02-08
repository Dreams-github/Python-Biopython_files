def orf(path):
    from Bio import SeqIO

    ss = SeqIO.parse(path, "fasta")

    for s in ss:
        f = open(s.id + ".fasta", "w")
        print(s.id)
        print("1st frame: ", s.seq.translate())
        print("2st frame: ", s.seq[1:].translate())
        print("3st frame: ", s.seq[2:].translate())
        rev = s.seq.reverse_complement()
        print("4th frame: ", rev.translate())
        print("5th frame: ", rev[1:].translate())
        print("6th frame: ", rev[2:].translate())
        out = ">" + s.id + "\n" + "f1:" + s.seq.translate() + "\n" + "f2:" + s.seq[1:].translate() + "\n" + "f3:" + s.seq[2:].translate() + "\n" + "f4:" + rev.translate() + "\n" + "f5:" + rev[
                                                                                                                                                                                1:].translate() + "\n" + "f6:" + rev[
                                                                                                                                                                                                                 2:].translate()
        print(out, file=f)
