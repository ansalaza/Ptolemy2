package gene_graph

import java.io.{File, PrintWriter}

import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile}
import utilities.MetaDataUtils.{Description, loadGeneFamilies}
import org.mapdb.{DBMaker, Serializer}
import utilities.SequenceGraphUtils.{empty_sequence_graph, empty_sequence_paths, msa2SequenceGraph, computeNodeEdgeCoverage}
import utilities.MSAutils.muscleMSA
import atk.ProgressBar.progress
import utilities.GFAutils.SequenceGraphWriter

/**
  * Author: Alex N. Salazar
  * Created on 10-3-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SequenceGraph extends SequenceGraphWriter {

  case class Config(
                     db: File = null,
                     geneFamilies: File = null,
                     outputDir: File = null,
                     prefix: String = null,
                     log: Int = 25,
                     path: String = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("sequence-graph") {
      opt[File]('s', "sequence-db") required() action { (x, c) =>
        c.copy(db = x)
      } text ("Sequence database as constructed by Ptolemy2.")
      opt[File]('f', "gene-families") required() action { (x, c) =>
        c.copy(geneFamilies = x)
      } text ("Gene-families file as constructed by Ptolemy2.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix name for output file.")
      note("\nOPTIONAL")
      opt[Int]("log-size") action { (x, c) =>
        c.copy(log = x)
      } text ("Report MSA progress after processing given number of gene-families (default is 25).")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.db)
      verifyFile(config.geneFamilies)
      sequenceGraph(config)
    }
  }

  def sequenceGraph(config: Config): Unit = {
    //load gene-families
    val gene_families = loadGeneFamilies(config.geneFamilies)
    println(timeStamp + "Found " + gene_families.size + " gene-family definitions")
    //load seq db
    val seqdb = DBMaker.fileDB(config.db).fileMmapEnable().make()
    //load hash map
    val seqmap = seqdb.treeMap("seqdb").keySerializer(Serializer.INTEGER).valueSerializer(Serializer.STRING).open()
    println(timeStamp + "Constructing sequence-graph (n=" + gene_families.size + "):")
    //iterate through gene families and construct sequence graph and paths
    val (graph, paths) = gene_families.foldLeft((empty_sequence_graph, empty_sequence_paths)) {
      case ((graph, paths), (family, genes)) => {
      progress(config.log)
      //set output file
      val local_seqs = new File(config.outputDir + "/" + family + ".fa")
      //set output fasta file
      val pw = new PrintWriter(local_seqs)
      //iterate through each gene sequence and output to file
      genes.foreach(gene => {
        pw.println(">" + gene)
        pw.println(seqmap.get(gene))
      })
      //close file
      pw.close()
      //run MSA
      val msa_file = muscleMSA(local_seqs, config.outputDir)
      //convert to graph
      val (local_graph, local_paths) = msa2SequenceGraph(msa_file, family)
      msa_file.delete()
      local_seqs.delete()
      //update global sequence graph, and global sequence paths
      (graph ++ local_graph, paths ++ local_paths)
    }}
    //close db
    seqmap.close()
    seqdb.close()
    println(timeStamp + "Computing coverages")
    //calculate node and edge coverages
    val (ncov, ecov) = computeNodeEdgeCoverage(graph, paths)
    println(timeStamp + "Writing to disk")
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".gfa")
    seqGraph2GFA(pw, graph, paths, ncov, ecov)
    pw.close()
    println(timeStamp + "Successfully completed")
  }

}
