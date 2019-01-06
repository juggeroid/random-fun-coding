object pwgen-scala { 
  def pwd(s: String = "", l: Int = 24) = {
           val rnd = new scala.util.Random((System.currentTimeMillis() + s.hashCode & ~(new scala.util.Random().nextInt() << 0xFFFF)).toLong)
           val str = (('a'to 'z').mkString) + ('A' to 'Z').mkString + ('0' to '9').mkString + "!@#$%^&*()-_=+/|{}[]~№;:?,<>"
           val r = Stream.continually((for (n <- 1 to l; c = str(rnd.nextInt(str.length))) yield c).mkString)
           r.head }; def main(args: Array[String]) = { for (n <- 1 to 24) println(pwd(n.toString, 24))}
}
