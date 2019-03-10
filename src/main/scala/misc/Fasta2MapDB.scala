package misc

import java.io.File

import org.mapdb.{DBMaker, Serializer}
import utilities.FileHandling.{getFileName, openFileWithIterator, verifyDirectory, verifyFile, timeStamp}
import utilities.SequenceFormatUtils.parseSeqName
import scala.annotation.tailrec
import utilities.MetaDataUtils.loadDescriptions

/**
  * Author: Alex N. Salazar
  * Created on 10-3-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object Fasta2MapDB {

  case class Config(
                     fasta: File = null,
                     outputDir: File = null,
                     descriptions: File = null
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("seq-db") {
      opt[File]('s', "sequences") required() action { (x, c) =>
        c.copy(fasta = x)
      } text ("All DNA/protein sequences in FASTA-format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      note("\nOPTIONAL")
      opt[File]("local-ids") action { (x, c) =>
        c.copy(descriptions = x)
      } text ("Sequence identifiers use local gene/protein IDs. Provide descriptions file to obtain mapping for " +
        "Ptolemy2-assigned IDs.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.fasta)
      seqDB(config)
    }
  }

  def seqDB(config: Config): Unit = {
    //load local-id -> gene ID
    val local2ID = {
      if(config.descriptions == null) Map[String, Int]()
      else loadDescriptions(config.descriptions).map(x => x.geneID -> x.id).toMap
    }
    //load file as iterator
    val iterator = openFileWithIterator(config.fasta)
    //load first fasta entry
    //get first sequence name
    val first_seqname = {
      assert(iterator.hasNext, "Expected to load first FASTA entry but found empty file")
      parseSeqName(iterator.next())
    }
    //create fileDB
    val db = DBMaker.fileDB(config.outputDir + "/" + getFileName(config.fasta) + ".db").make()
    //create map
    var map = db.treeMap("seqdb").keySerializer(Serializer.INTEGER).valueSerializer(Serializer.STRING).createOrOpen()

    @tailrec def updateSeqDB(name: String, seq: StringBuilder): Unit = {
      def updateMap(): Unit = {
        val id: Int = local2ID.getOrElse(name, name.toInt)
        map.put(id, seq.toString)
      }
      //iterator is empty, attempt to fetch last entry
      if (iterator.isEmpty) {
        assert(seq.nonEmpty, "Expected sequence for last FASTA entry:" + name)
        updateMap
      }
      else {
        //end of file
        if (!iterator.hasNext) {
          //sanity check
          assert(seq.nonEmpty, "Expected sequence on last FASTA entry")
          //update
          updateMap
        }
        else {
          //get next line
          val next = iterator.next
          //continuing current FASTA entry
          if (!next.startsWith(">")) updateSeqDB(name, seq.append(next))
          //start of new FASTA entry
          else {
            //update
            updateMap
            //move on
            updateSeqDB(parseSeqName(next), new StringBuilder)
          }
        }
      }
    }
    println(timeStamp + "Processing sequences and creating database")
    //iterate through fasta file and update map
    updateSeqDB(first_seqname, new StringBuilder)
    //close map
    map.close
    println(timeStamp + "Sucessfully completed!")
  }

}
