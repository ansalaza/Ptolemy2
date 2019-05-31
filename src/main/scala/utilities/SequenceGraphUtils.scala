package utilities

import java.io.File

import utilities.FileHandling.{openFileWithIterator, timeStamp}
import utilities.SequenceFormatUtils.parseSeqName
import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 10-3-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SequenceGraphUtils {

  case class AminoAcid(seq_id: Int, gene_id: Int, aa: Char) {
    def gfaID(): String = gene_id + "_" + seq_id + "_" + aa
  }

  type SequenceGraph = Map[AminoAcid, List[AminoAcid]]
  type SequencePaths = Map[String, List[AminoAcid]]
  val empty_sequence_graph: SequenceGraph = Map[AminoAcid, List[AminoAcid]]()
  val empty_sequence_paths: SequencePaths = Map[String, List[AminoAcid]]()


  /**
    * Function to parse through FASTA-formatted MSA file and output SequenceGraph and SequencePaths objects
    *
    * @return 2-tuple: (SequenceGraph, SequencePaths)
    */
  def msa2SequenceGraph(msa_file: File, gene_id: Int): (SequenceGraph, SequencePaths) = {
    //load msa file as iterator
    val iterator = openFileWithIterator(msa_file)
    //get first sequence name
    val first_seqname = {
      assert(iterator.hasNext, "Expected to load first FASTA entry but found empty file")
      parseSeqName(iterator.next())
    }


    /**
      * Tail-recursive method to iterate through MSA FASTA file and create a sequence graph and paths.
      *
      * @param seqname Current sequende name
      * @param seq     Accumulating sequence
      * @param graph   Sequence graph
      * @param paths   Sequence path
      * @return 2-tuple: (SequenceGraph, SequencePaths)
      */
    @tailrec def _msa2SequenceGraph(seqname: String,
                                    seq: StringBuilder,
                                    graph: SequenceGraph,
                                    paths: SequencePaths): (SequenceGraph, SequencePaths) = {

      /**
        * Function to create an amino acid sequence using the current seq
        *
        * @return List[AminoAcid]
        */
      def setSequence(): List[AminoAcid] =
        seq.toString().toList.zipWithIndex.filter(_._1 != '-').map(x => new AminoAcid(x._2, gene_id, x._1))

      /**
        * Function to update sequence graph given a sequence
        *
        * @param sequence Sequence
        * @return Updated SequenceGraph
        */
      def updateGraph(sequence: List[AminoAcid]): SequenceGraph = {
        //update according to sequence size
        sequence.size match {
          case 1 => {
            println(timeStamp + "Found sequence of one amino acid: " + (seqname, sequence))
            //update graph
            graph + (sequence.head -> graph.getOrElse(sequence.head, List[AminoAcid]()))
          }
          case _ => {
            //iterate through sequence and update graph
            sequence.sliding(2).map(x => (x(0), x(1))).foldLeft(graph)((acc, edge) => {
              acc + (edge._1 -> (edge._2 :: acc.getOrElse(edge._1, List[AminoAcid]()))) +
                (edge._2 -> (acc.getOrElse(edge._2, List[AminoAcid]())))
            })
          }
        }
      }

      /**
        * Function to update path sequence
        *
        * @param sequence Sequence
        * @return UpdatedPath
        */
      def updatePath(sequence: List[AminoAcid]): SequencePaths = paths + (seqname -> sequence)

      //end of file
      if (!iterator.hasNext) {
        //sanity check
        assert(seq.nonEmpty, "Expected sequence on last FASTA entry")
        //set sequence
        val sequence = setSequence()
        //update
        (updateGraph(sequence).mapValues(_.distinct), updatePath(sequence))
      }
      else {
        //get next line
        val next = iterator.next
        //continuing current FASTA entry
        if (!next.startsWith(">")) _msa2SequenceGraph(seqname, seq.append(next), graph, paths)
        //start of new FASTA entry
        else {
          //set sequence
          val sequence = setSequence()
          //update and move on
          _msa2SequenceGraph(parseSeqName(next), new StringBuilder, updateGraph(sequence), updatePath(sequence))
        }
      }
    }

    _msa2SequenceGraph(first_seqname, new StringBuilder, Map(), Map())
  }


  def computeNodeEdgeCoverage(graph: SequenceGraph,
                              paths: SequencePaths
                             ): (Map[String, Int], Map[(String, String), Int]) = {
    //compute node coverage
    val ncov = graph.keys.toList.groupBy(_.gfaID()).mapValues(_.size)
    //compute edge coverage
    val ecov = paths.toList.foldLeft(Map[(String, String), Int]()) { case (cov, (name, path)) => {
      path.size match {
        //no edges, move on
        case 1 => cov
        //process every forward edge
        case _ => path.sliding(2).map(x => (x(0), x(1))).foldLeft(cov)((acc, e) =>{
          assert(graph.contains(e._1), "Expected node " + e._1 + " from path " + name + " to exist in graph: " + path)
          assert(graph.contains(e._2), "Expected node " + e._2 + " from path " + name + " to exist in graph: " + path)
          acc + ((e._1.gfaID(), e._2.gfaID()) -> (1 + acc.getOrElse((e._1.gfaID(), e._2.gfaID()), 0)))
        })
      }
    }
    }
    (ncov, ecov)
  }

}
