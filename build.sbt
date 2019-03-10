name := "Ptolemy2"

version := "0.1"

scalaVersion := "2.12.8"

libraryDependencies ++= Seq(
  "org.scalatest" %% "scalatest" % "3.0.4" % "test",
  "com.github.scopt" % "scopt_2.12" % "3.7.0",
  "org.mapdb" % "mapdb" % "3.0.7",
  "de.sciss" %% "fingertree" % "1.5.2"
)