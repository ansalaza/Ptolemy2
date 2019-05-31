package utilities

import utilities.FileHandling.timeStamp

/**
  * Author: Alex N. Salazar
  * Created on 4-4-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SequenceUtils {

  /**
    * Function to convert a given DNA sequence into protein sequence
    *
    * @return String
    */
  def dna2Protein(seq: String, showWarnings: Boolean = false): String = {
    val codons = seq.grouped(3).toList.zipWithIndex
    val ct = codons.size - 1
    //codons.foreach(x => assert(x.size == 3, "Unexpected codon size " + x.size + " in " +  codons))
    codons.map(codon => {
      //try to fetch protein
      val prot = codon2protein.get(codon._1)
      if (prot.isEmpty || prot.get == '_') {
        if(showWarnings){
          if (prot.isEmpty) println(timeStamp + "Unexpected codon sequence " + codon + ". " + codons)
          if (prot.nonEmpty && prot.get == '_' && codon._2 != ct)
            println(timeStamp + "Unexpected stopping codon " + codon + ". " + codons)
        }
        ""
      } else prot.get
    }).mkString("")
  }

  //TODO: account for sequence with ambiguous/non-ACTG bases
  /**
    * Function to obtain reverse complement for a given DNA sequence
    *
    * @return String
    */
  def reverseComplement: String => String = seq => seq.reverse.map(nt => {
    nt match {
      case 'A' => 'T'
      case 'T' => 'A'
      case 'C' => 'G'
      case 'G' => 'C'
      case _ => 'N'
    }
  })

  private val codon2protein = Map(
    "ATA" -> 'I', "ATC" -> 'I', "ATT" -> 'I', "ATG" -> 'M', "ACA" -> 'T', "ACC" -> 'T', "ACG" -> 'T', "ACT" -> 'T',
    "AAC" -> 'N', "AAT" -> 'N', "AAA" -> 'K', "AAG" -> 'K', "AGC" -> 'S', "AGT" -> 'S', "AGA" -> 'R', "AGG" -> 'R', "CTA" -> 'L', "CTC" -> 'L',
    "CTG" -> 'L', "CTT" -> 'L', "CCA" -> 'P', "CCC" -> 'P', "CCG" -> 'P', "CCT" -> 'P', "CAC" -> 'H', "CAT" -> 'H', "CAA" -> 'Q', "CAG" -> 'Q',
    "CGA" -> 'R', "CGC" -> 'R', "CGG" -> 'R', "CGT" -> 'R', "GTA" -> 'V', "GTC" -> 'V', "GTG" -> 'V', "GTT" -> 'V', "GCA" -> 'A', "GCC" -> 'A',
    "GCG" -> 'A', "GCT" -> 'A', "GAC" -> 'D', "GAT" -> 'D', "GAA" -> 'E', "GAG" -> 'E', "GGA" -> 'G', "GGC" -> 'G', "GGG" -> 'G', "GGT" -> 'G',
    "TCA" -> 'S', "TCC" -> 'S', "TCG" -> 'S', "TCT" -> 'S', "TTC" -> 'F', "TTT" -> 'F', "TTA" -> 'L', "TTG" -> 'L', "TAC" -> 'Y', "TAT" -> 'Y',
    "TAA" -> '_', "TAG" -> '_', "TGC" -> 'C', "TGT" -> 'C', "TGA" -> '_', "TGG" -> 'W')

}
