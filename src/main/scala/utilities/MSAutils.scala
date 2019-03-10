package utilities

import java.io.File
import utilities.FileHandling.getFileName
import scala.sys.process.{Process, ProcessLogger}

/**
  * Author: Alex N. Salazar
  * Created on 10-3-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object MSAutils {

  def muscleMSA(input: File, outputdir: File): File = {
    /*
     Capture stdout and stderr.
     Note: using mutable stringbuilder, dirty but gets job done; may be optimizedlater
    */
    var out = new StringBuilder
    var err = new StringBuilder
    val logger = ProcessLogger((o: String) => out.append(o + "\n"), (e: String) => err.append(e + "\n"))
    //set output file
    val output_file = new File(outputdir + "/" + getFileName(input) + ".msa.fa")
    //set command
    val command = Seq("muscle", "-in", input.getAbsolutePath, "-out", output_file.getAbsolutePath)
    //run alignment command
    val exitcode = Process(command).!(logger)
    //sanity check
    assert(exitcode == 0, "Non-zero exitcode: " + err.toString())
    output_file
  }

}
