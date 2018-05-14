package simulations

/**
 * This class defines the information for various "cuts" to the spheres.
 */
trait CollisionForce extends ((MutableBody, MutableBody) => Unit)
trait OptimizedCollisionForce extends ((MutableBody, MutableBody, MVect3, Double) => Unit)

class SphereCollide(b1: Double, damping: Double) extends CollisionForce {
  val opt = new OptimizedSphereCollide(b1, damping)
  def apply(pi: MutableBody, pj: MutableBody): Unit = {
    val dx = pi.p.x - pj.p.x
    val dy = pi.p.y - pj.p.y
    val dz = pi.p.z - pj.p.z
    val distSqr = dx * dx + dy * dy + dz * dz
    if (distSqr < (pi.radius + pj.radius) * (pi.radius + pj.radius)) {
      opt(pi, pj, new MVect3(dx, dy, dz), math.sqrt(distSqr))
    }
  }
}

class OptimizedSphereCollide(b1: Double, damping: Double) extends OptimizedCollisionForce {
  def apply(pi: MutableBody, pj: MutableBody, n: MVect3, dist: Double): Unit = {
    val overlap = (pi.radius + pj.radius) - dist
    n.normalize
    val v = -(n dot (pi.v - pj.v))
    val mag = -b1 * overlap - v * damping
    pi.a.x -= n.x * mag // pi.mass
    pi.a.y -= n.y * mag // pi.mass
    pi.a.z -= n.z * mag // pi.mass
    pj.a.x += n.x * mag // pj.mass
    pj.a.y += n.y * mag // pj.mass
    pj.a.z += n.z * mag // pj.mass
  }
}

case class SofteningData(sd: MVect3, cutoff: Double, soft: Double)

/**
 * This is the forcing for a collision where the particle's restoring force has been softened in a particular direction.
 */
class SoftenedSphereCollide(b1: Double, damping: Double, soft: Seq[SofteningData]) extends CollisionForce {
  val opt = new OptimizedSoftenedSphereCollide(b1, damping, soft)
  def apply(pi: MutableBody, pj: MutableBody): Unit = {
    val dx = pi.p.x - pj.p.x
    val dy = pi.p.y - pj.p.y
    val dz = pi.p.z - pj.p.z
    val distSqr = dx * dx + dy * dy + dz * dz
    if (distSqr < (pi.radius + pj.radius) * (pi.radius + pj.radius)) {
      opt(pi, pj, new MVect3(dx, dy, dz), math.sqrt(distSqr))
    }
  }
}

class OptimizedSoftenedSphereCollide(b1: Double, damping: Double, soft: Seq[SofteningData]) extends OptimizedCollisionForce {
  def f(x: Double, cutoff: Double): Double = (x - cutoff) / (1 - cutoff) // fixxxxxxxxxxx

  def apply(pi: MutableBody, pj: MutableBody, n: MVect3, dist: Double): Unit = {
    val overlap = (pi.radius + pj.radius) - dist
    n.normalize
    val closestSofti = soft.maxBy(_.sd dot n)
    val closestSoftj = soft.minBy(_.sd dot n)
    val softFactori = closestSofti.sd dot n // This need to be changed to handle multiple directions and pick the right one.
    val softFactorj = closestSoftj.sd dot n // This need to be changed to handle multiple directions and pick the right one.
    val softBi = if (softFactori > closestSofti.cutoff) b1 * (1 - closestSofti.soft * f(softFactori, closestSofti.cutoff)) else b1
    val softBj = if (softFactorj > closestSoftj.cutoff) b1 * (1 - closestSoftj.soft * f(softFactorj, closestSoftj.cutoff)) else b1

    val v = -(n dot (pi.v - pj.v))
    val mag = -softBi * softBj * overlap - v * damping
    pi.a.x -= n.x * mag // pi.mass
    pi.a.y -= n.y * mag // pi.mass
    pi.a.z -= n.z * mag // pi.mass
    pj.a.x += n.x * mag // pj.mass
    pj.a.y += n.y * mag // pj.mass
    pj.a.z += n.z * mag // pj.mass
  }
}

case class WarpingData(sd: MVect3, cutoff: Double, dentFactor: Double)

/**
 * This is the forcing for a collision where the particle's radius is modified in a particular direction.
 */
class WarpedCollide(b1: Double, damping: Double, warp: Seq[WarpingData]) extends CollisionForce {
  val opt = new OptimizedWarpedCollide(b1, damping, warp)
  def apply(pi: MutableBody, pj: MutableBody): Unit = {
    val dx = pi.p.x - pj.p.x
    val dy = pi.p.y - pj.p.y
    val dz = pi.p.z - pj.p.z
    val distSqr = dx * dx + dy * dy + dz * dz
    if (distSqr < (pi.radius + pj.radius) * (pi.radius + pj.radius)) {
      opt(pi, pj, new MVect3(dx, dy, dz), math.sqrt(distSqr))
    }
  }
}

class OptimizedWarpedCollide(b1: Double, damping: Double, warp: Seq[WarpingData]) extends OptimizedCollisionForce {
  def apply(pi: MutableBody, pj: MutableBody, n: MVect3, dist: Double): Unit = {
    n.normalize
    val closestWarpi = warp.maxBy(_.sd dot n)
    val closestWarpj = warp.minBy(_.sd dot n)
    val iFactor = closestWarpi.sd dot n // These need to be changed to handle multiple directions and pick the right one.
    val jFactor = -(closestWarpj.sd dot n)
    val pi0 = pi.radius
    val pi1 = pi0 * closestWarpi.dentFactor
    val cutoffi = closestWarpi.cutoff
    val irad = if (iFactor > cutoffi) pi0 + (iFactor - cutoffi) / (1.0 - cutoffi) * (pi1 - pi0) else pi0
    val pj0 = pj.radius
    val pj1 = pj0 * closestWarpj.dentFactor
    val cutoffj = closestWarpj.cutoff
    val jrad = if (jFactor > cutoffj) pj0 + (jFactor - cutoffj) / (1.0 - cutoffj) * (pj1 - pj0) else pj0

    if (dist < irad + jrad) {
      val overlap = (irad + jrad) - dist

      val v = -(n dot (pi.v - pj.v))
      val mag = -b1 * overlap - v * damping
      pi.a.x -= n.x * mag // pi.mass
      pi.a.y -= n.y * mag // pi.mass
      pi.a.z -= n.z * mag // pi.mass
      pj.a.x += n.x * mag // pj.mass
      pj.a.y += n.y * mag // pj.mass
      pj.a.z += n.z * mag // pj.mass
    }
  }
}
