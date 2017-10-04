package simulations

/**
 * The purpose of this code is to test the coefficient of restitution and related values for simple systems where
 * two bodies collide with one another. We want to explore a high dimensional parameter space that includes:
 *  - b1: constant determining particle hardness
 *  - damping: constant determining the damping of particle velocities during overlap
 *  - dt: time step
 *  - impact velocity
 *  
 * For the non-spherical collisions we also have the cutoff and the softness/dent factors.
 */
object CoefRestTesting extends App {
  
  
  def runTest(p1: MutableBody, p2: MutableBody, dt: Double): (Double, Double, Double) = {
    ???
  }
}