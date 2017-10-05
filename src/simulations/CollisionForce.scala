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

/**
 * This is the forcing for a collision where the particle's restoring force has been softened in a particular direction.
 */
class SoftenedSphereCollide(b1: Double, damping: Double, sd: MVect3, cutoff: Double, soft: Double) extends CollisionForce {
  val opt = new OptimizedSoftenedSphereCollide(b1, damping, sd, cutoff, soft)
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

class OptimizedSoftenedSphereCollide(b1: Double, damping: Double, sd: MVect3, cutoff: Double, soft: Double) extends OptimizedCollisionForce {
  def f(x: Double): Double = (x - cutoff) / (1 - cutoff) // fixxxxxxxxxxx

  def apply(pi: MutableBody, pj: MutableBody, n: MVect3, dist: Double): Unit = {
    val overlap = (pi.radius + pj.radius) - dist
    n.normalize
    val softFactor = sd dot n // This need to be changed to handle multiple directions and pick the right one.
    val softB = if (softFactor > cutoff) b1 * (1 - soft * f(softFactor)) else b1

    val v = -(n dot (pi.v - pj.v))
    val mag = -softB * overlap - v * damping
    pi.a.x -= n.x * mag // pi.mass
    pi.a.y -= n.y * mag // pi.mass
    pi.a.z -= n.z * mag // pi.mass
    pj.a.x += n.x * mag // pj.mass
    pj.a.y += n.y * mag // pj.mass
    pj.a.z += n.z * mag // pj.mass
  }
}

/**
 * This is the forcing for a collision where the particle's radius is modified in a particular direction.
 */
class WarpedCollide(b1: Double, damping: Double, sd: MVect3, cutoff: Double, dentFactor: Double) extends CollisionForce {
  val opt = new OptimizedWarpedCollide(b1, damping, sd, cutoff, dentFactor)
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

class OptimizedWarpedCollide(b1: Double, damping: Double, sd: MVect3, cutoff: Double, dentFactor: Double) extends OptimizedCollisionForce {
  def apply(pi: MutableBody, pj: MutableBody, n: MVect3, dist: Double): Unit = {
    n.normalize
    val iFactor = sd dot n // These need to be changed to handle multiple directions and pick the right one.
    val jFactor = -(sd dot n)
    val pi0 = pi.radius
    val pi1 = pi0 * dentFactor
    val irad = if (iFactor > cutoff) pi0 + (iFactor - cutoff) / (1.0 - cutoff) * (pi1 - pi0) else pi0
    val pj0 = pj.radius
    val pj1 = pj0 * dentFactor
    val jrad = if (jFactor > cutoff) pj0 + (jFactor - cutoff) / (1.0 - cutoff) * (pj1 - pj0) else pj0

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
