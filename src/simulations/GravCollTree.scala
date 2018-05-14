package simulations

import collection.GenSeq

/**
 * I'm going to use a quadtree here because I think that will work well for this particular application.
 * Things won't get too unbalanced, and I am going to stop splitting at the size of the largest particle,
 * so the tree has a fixed maximum depth of ~log2(cellSize/particleSize).
 */
class GravCollTree(p: Array[MutableBody], size: Double) {
  import GravCollTree._

  val maxRad = p.map(_.radius).max
  val root = buildTree(p, 0.0, 0.0, size)

  val theta = 0.5

  def addAccel(pi: MutableBody, collForcing: (MutableBody, MutableBody, MVect3, Double) => Unit, offx: Double, offy: Double): Unit = {
    def helper(n: Node): Unit = {
      val dx = pi.p.x - (n.cx+offx)
      val dy = pi.p.y - (n.cy+offy)
      val dist = math sqrt (dx * dx + dy * dy)
      if (dist * theta > n.size) {
        singleAccelCM(pi, n.cmData.x+offx, n.cmData.y+offy, n.cmData.z, n.cmData.mass)
      } else if (n.bodies.nonEmpty) {
        for (o <- n.bodies) if (o != pi) {
          singleAccelPart(pi, o, collForcing, offx, offy)
        }
      } else {
        n.children.foreach(helper)
      }
    }
    helper(root)
  }

  private def singleAccelCM(pi: MutableBody, ox: Double, oy: Double, oz: Double, omass: Double): Unit = {
    val dx = pi.p.x - ox
    val dy = pi.p.y - oy
    val dz = pi.p.z - oz
    val dist = math.sqrt(dx * dx + dy * dy + dz * dz)
    val magi = omass / (dist * dist * dist)
    pi.a.x -= magi * dx
    pi.a.y -= magi * dy
    pi.a.z -= magi * dz
  }

  private def singleAccelPart(pi: MutableBody, pj: MutableBody, collForcing: (MutableBody, MutableBody, MVect3, Double) => Unit,
      offx: Double, offy: Double): Unit = {
    val dx = pi.p.x - (pj.p.x+offx)
    val dy = pi.p.y - (pj.p.y+offy)
    val dz = pi.p.z - pj.p.z
    val dist = math.sqrt(dx * dx + dy * dy + dz * dz)
    val magi = pj.mass / (dist * dist * dist)
    pi.a.x -= magi * dx
    pi.a.y -= magi * dy
    pi.a.z -= magi * dz
    if (dist < pi.radius+pj.radius && pi.hashCode() < pj.hashCode()) 
      collForcing(pi, pj, new MVect3(dx, dy, dz), dist)
  }

  private def buildTree(parts: GenSeq[MutableBody], cx: Double, cy: Double, size: Double): Node = {
    if (size / 2 < maxRad) Node(emptyChildren, cx, cy, size, parts.toArray, calcCMsP(parts)) else {
      val groups = parts.par.groupBy(childIndex(_, cx, cy))
      val s2 = size / 2
      val s4 = size / 4
      val children = (0 to 3).par.map { c =>
        groups.get(c).map { ps =>
          val ccx = cx + (if ((c & 1) == 0) -s4 else s4)
          val ccy = cy + (if ((c & 2) == 0) -s4 else s4)
          buildTree(ps, ccx, ccy, s2)
        }.getOrElse(Node(emptyChildren, cx, cy, size, emptyBodies, emptyCM))
      }.toArray
      Node(children, cx, cy, size, emptyBodies, calcCMs(children))
    }
  }

  private def childIndex(pi: MutableBody, cx: Double, cy: Double): Int = {
    (if (pi.p.x > cx) 1 else 0) + (if (pi.p.y > cy) 2 else 0)
  }

  private def calcCMsP(parts: GenSeq[MutableBody]): CMData = {
    val (xsum, ysum, zsum, m) = parts.foldLeft((0.0, 0.0, 0.0, 0.0)) {
      case ((xs, ys, zs, ms), ps) =>
        (xs + ps.p.x * ps.mass, ys + ps.p.y * ps.mass, zs + ps.p.z * ps.mass, ms + ps.mass)
    }
    CMData(xsum / m, ysum / m, zsum / m, m)
  }

  private def calcCMs(children: GenSeq[Node]): CMData = {
    val (xsum, ysum, zsum, m) = children.foldLeft((0.0, 0.0, 0.0, 0.0)) {
      case (t, null) => t 
      case ((xs, ys, zs, ms), n) =>
        (xs + n.cmData.x * n.cmData.mass, ys + n.cmData.y * n.cmData.mass, zs + n.cmData.z * n.cmData.mass, ms + n.cmData.mass)
    }
    CMData(xsum / m, ysum / m, zsum / m, m)
  }
}

object GravCollTree {
  val emptyChildren = Array.empty[Node]
  val emptyBodies = Array.empty[MutableBody]
  val emptyCM = CMData(0.0, 0.0, 0.0, 0.0)
  case class CMData(x: Double, y: Double, z: Double, mass: Double)
  case class Node(children: Array[Node], cx: Double, cy: Double, size: Double, bodies: Array[MutableBody],
                  cmData: CMData)
}
