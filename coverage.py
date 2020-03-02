import unittest

reads_path = "reads.csv"
loci_path = "loci.csv"
lowest_multiple_of_ten = 10000

def load_reads():
    """Returns a dictionary containing the reads.
       Format: {multiple of lowest_multiple of ten: [(start, length),...]}"""

    reads = {}
    with open(reads_path, "r", encoding="utf-8") as fo:
        next(fo) # skip headers
        for line in fo:
            line = line.split(",")
            start = int(line[0])
            length = int(line[1])
            key = start // lowest_multiple_of_ten
            read_tuple = (start, length)

            if reads.get(key) is None:      # creates a new key-value pair if needed
                reads[key] = [read_tuple]
            else:                           # else appends to existing value
                reads[key].append(read_tuple)

            # a read will also need to be appended to the lists of higher multiples
            # if they are covered by the range its length goes over
            key += 1
            while start + length > key*lowest_multiple_of_ten:
                if reads.get(key) is None:      # creates a new key-value pair if needed
                    reads[key] = [read_tuple]
                else:                           # else appends to existing value
                    reads[key].append(read_tuple)
                key += 1
    return reads

def load_loci():
    """Returns a list containing the loci."""

    loci = []
    with open(loci_path, "r", encoding="utf-8") as fo:
        next(fo) # skip headers
        for line in fo:
            loci.append(int(line.split(",")[0]))
    return loci

def save_coverages(loci, coverages):
    """Saves coverages to loci_path, given a list of loci and a list of coverages."""

    with open(loci_path, "w", encoding="utf-8") as fo:
        fo.write("position,coverage\n")
        fo.writelines([f"{locus},{coverage}\n" for locus, coverage in zip(loci, coverages)])

def is_covered(read, position):
    """Returns true if position is covered by read, otherwise false."""

    return position >= read[0] and position < read[0]+read[1]

def calculate_coverages(reads, loci):
    """Returns a list of the coverage for each position in loci."""

    coverages = []
    for locus in loci:
        count = 0
        curr_reads = reads.get(locus // lowest_multiple_of_ten)
        if curr_reads is not None:
            for read in curr_reads:
                if is_covered(read, locus):
                    count += 1
        coverages.append(count)
    return coverages

class TestCoverage(unittest.TestCase):

    def test_load_reads(self):
        reads = load_reads()
        self.assertEqual(len(reads[5]), 2)
        self.assertIn((50532, 125), reads[5])

    def test_calculate_coverages(self):
        expected_output = [1190, 12206, 141295, 6809, 6993, 1096, 62732,
                           43740, 78510, 42854, 62, 221, 24560, 1, 2581,
                           84, 6209, 36, 431, 213040]
        reads = load_reads()
        loci = load_loci()
        calculated_output = calculate_coverages(reads, loci)

        for i in range(len(expected_output)):
            self.assertEqual(expected_output[i], calculated_output[i])

if __name__ == "__main__":
    # unittest.main()   # uncomment to test first 20 loci

    reads = load_reads()
    loci = load_loci()

    save_coverages(loci, calculate_coverages(reads, loci))
