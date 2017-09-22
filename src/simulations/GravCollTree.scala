package simulations

/**
 * I'm going to use a quadtree here because I think that will work well for this particular application.
 * Things won't get too unbalanced, and I am going to stop splitting at the size of the largest particle,
 * so the tree has a fixed maximum depth of ~log2(cellSize/particleSize).
 */
class GravCollTree(p: Array[MutableBody], size: Double) {
	import GravCollTree._

	val maxRad = p.map(_.radius).max
  val root = buildTree(p, 0.0, 0.0, size)
  
  val theta = 0.3

  def addAccel(pi: MutableBody, collForcing: (MutableBody, MutableBody) => Unit): Unit = {
    def helper(n: Node): Unit = {
//      println("helper on "+n)
      val dx = pi.p.x - n.cx
      val dy = pi.p.y - n.cy
      val dist = math sqrt(dx*dx + dy*dy)
      if(dist*theta > n.size) {
        singleAccel(pi, n.cmData.x, n.cmData.y, n.cmData.z, n.cmData.mass)
      } else if(n.bodies.nonEmpty) {
        for(o <- n.bodies; if o != pi) {
          singleAccel(pi, o.p.x, o.p.y, o.p.z, o.mass)
          if(pi.hashCode() < o.hashCode()) collForcing(pi, o)
        }
      } else {
        n.children.foreach(helper)
      }
    }
    helper(root)
  }
  
  private def singleAccel(pi: MutableBody, ox: Double, oy: Double, oz: Double, omass: Double): Unit = {
    val dx = pi.p.x - ox
    val dy = pi.p.y - oy
    val dz = pi.p.z - oz
    val dist = math.sqrt(dx * dx + dy * dy + dz * dz)
    val magi = omass / (dist * dist * dist)
    pi.a.x -= magi * dx
    pi.a.y -= magi * dy
    pi.a.z -= magi * dz
  }

  private def buildTree(parts: Array[MutableBody], cx: Double, cy: Double, size: Double): Node = {
    if (size / 2 < maxRad) Node(emptyChildren, cx, cy, size, parts, calcCMs(parts)) else {
      val groups = parts.groupBy(childIndex(_, cx, cy))
      val s2 = size / 2
      val s4 = size / 4
      val children = groups.toArray.sortBy(_._1).map {
        case (c, ps) =>
          val ccx = cx + (if ((c & 1) == 0) -s4 else s4)
          val ccy = cy + (if ((c & 2) == 0) -s4 else s4)
          buildTree(ps, ccx, ccy, s2)
      }
      Node(children, cx, cy, size, emptyBodies, calcCMs(children))
    }
  }

  private def childIndex(pi: MutableBody, cx: Double, cy: Double): Int = {
    (if (pi.p.x > cx) 1 else 0) + (if (pi.p.y > cy) 2 else 0)
  }

  private def calcCMs(parts: Array[MutableBody]): CMData = {
    val (xsum, ysum, zsum, m) = parts.foldLeft((0.0, 0.0, 0.0, 0.0)) {
      case ((xs, ys, zs, ms), ps) =>
        (xs + ps.p.x * ps.mass, ys + ps.p.y * ps.mass, zs + ps.p.z * ps.mass, ms + ps.mass)
    }
    CMData(xsum / m, ysum / m, zsum / m, m)
  }

  private def calcCMs(children: Array[Node]): CMData = {
    val (xsum, ysum, zsum, m) = children.foldLeft((0.0, 0.0, 0.0, 0.0)) {
      case ((xs, ys, zs, ms), n) =>
        (xs + n.cmData.x * n.cmData.mass, ys + n.cmData.y * n.cmData.mass, zs + n.cmData.z * n.cmData.mass, ms + n.cmData.mass)
    }
    CMData(xsum / m, ysum / m, zsum / m, m)
  }
}

object GravCollTree {
  val emptyChildren = Array.empty[Node]
  val emptyBodies = Array.empty[MutableBody]
  case class CMData(x: Double, y: Double, z: Double, mass: Double)
  case class Node(children: Array[Node], cx: Double, cy: Double, size: Double, bodies: Array[MutableBody],
                  cmData: CMData)
}