
module Main where

import Graphics.Gloss
import Graphics.Gloss.Data.ViewPort
import qualified Graphics.Gloss.Data.Point.Arithmetic as PA
import Graphics.Gloss.Data.Vector (magV, rotateV, argV, normalizeV, dotV)
import Graphics.Gloss.Geometry.Angle (radToDeg)

import System.Random
import Data.List (unfoldr)
import System.Console.ParseArgs

type Position = Vector
type Velocity = Vector
type Force = Vector

data Boid = Boid Position Velocity deriving (Eq)

polarToCartesian :: Float -> Float -> Vector
polarToCartesian r theta = (r * cos theta, r * sin theta)

randomRange :: (RandomGen g, UniformRange a) => a -> a -> g -> ([a], g)
randomRange lower upper g = (unfoldr (Just . uniformR (lower, upper)) g, g)

randomVectors :: RandomGen g => Float -> g -> ([Vector], g)
randomVectors r g =
    let (rs, g1) = randomRange 0.0 r g
        (alphas, g2) = randomRange 0.0 (2*pi :: Float) g1 in
            (zipWith polarToCartesian rs alphas, g2)

newtype Boids = Boids [Boid]

generateBoids :: RandomGen g => Int -> Int -> g -> Boids
generateBoids n spacing g =
    let points = [(fromIntegral $ x * spacing,
                   fromIntegral $ y * spacing)
                   | x <- [0..n-1], y <- [0..n-1]]
        (velocities, _) = randomVectors 100 g in
            Boids $ zipWith Boid points velocities

getCentralPoint :: Boids -> Position
getCentralPoint (Boids boids) =
    let boidPositions = map (\(Boid p' _) -> p') boids in
        (1 / fromIntegral (length boids)) PA.* foldr (PA.+) (0, 0) boidPositions

centeredBoids :: Boids -> Boids
centeredBoids bs@(Boids boids) = Boids . map (\(Boid p v) -> Boid (p PA.- c) v) $ boids
    where c = getCentralPoint bs

model :: RandomGen g => Int -> Int -> g -> Boids
model numBoidsPerSide boidInitialSpacing g = centeredBoids $ generateBoids numBoidsPerSide boidInitialSpacing g

boidVertices :: Path
boidVertices = [(-1, -2), (0, 2), (1, -2)]

drawBoid :: Boid -> Picture
drawBoid boid@(Boid p v) =
    uncurry Translate p . Scale 3 3 . Rotate (90 - (radToDeg . argV $ v))
        . Polygon $ boidVertices

drawGrid :: Int -> Picture
drawGrid spacing = Color (greyN 0.8) . Pictures $
    [Line [(x, -10000), (x, 10000)] | x <- map (fromIntegral . (* spacing)) [-1000..1000]] ++
    [Line [(-10000, y), (10000, y)] | y <- map (fromIntegral . (* spacing)) [-1000..1000]]

drawModel :: Boids -> Picture
drawModel bs@(Boids boids) =
    uncurry Translate (PA.negate c) . Pictures $ drawGrid 50 : map drawBoid boids
        where c = getCentralPoint bs

getNearbyBoids :: SimulationParams -> Boids -> Boid -> Boids
getNearbyBoids params (Boids boids) b1@(Boid p _) =
    Boids . filter (\b2@(Boid p2 _) ->
                        b1 /= b2 && magV (p2 PA.- p) < nearbyRadius params)
                $ boids

separation :: SimulationParams -> Boids -> Boid -> Force
separation params bs@(Boids boids) b@(Boid p v) =
    let Boids nearbyBoids = getNearbyBoids params bs b in
        foldr ((PA.+) .
                (\(Boid p2 _) ->
                    let d = p PA.- p2 in
                        (separationStrength params / dotV d d) PA.* normalizeV d))
            (0, 0) nearbyBoids

average :: Fractional a => [a] -> a
average b = sum b / fromIntegral (length b)

alignment :: SimulationParams -> Boids -> Boid -> Force
alignment params bs@(Boids boids) b@(Boid _ v) =
    let Boids nearbyBoids = getNearbyBoids params bs b in
        if null nearbyBoids then
            (0, 0)
        else
            let averageAlignment = average .
                                    map (\(Boid _ v') -> argV v') $ nearbyBoids
                currentAlignment = argV v
                difference = currentAlignment - averageAlignment in
                    signum difference * max (abs difference) (maxAlignmentSpeed params)
                        PA.* rotateV (-pi/2) v

cohesion :: SimulationParams -> Boids -> Boid -> Force
cohesion params bs@(Boids boids) b@(Boid p _) =
    let Boids nearbyBoids = getNearbyBoids params bs b in
        if null nearbyBoids then
            (0, 0)
        else
            let centreOfMass = getCentralPoint bs in
                cohesionStrength params PA.* (centreOfMass PA.- p)

calculateForce :: SimulationParams -> Boids -> Boid -> Force
calculateForce params boids boid =
    let force = separation params boids boid PA.+ alignment params boids boid PA.+ cohesion params boids boid in
        if magV force > maxForce params then
            maxForce params PA.* normalizeV force
        else
            force

stepModel :: SimulationParams -> ViewPort -> Float -> Boids -> Boids
stepModel params v f bs@(Boids boids) =
    Boids . map (\b@(Boid p v) ->
        Boid (p PA.+ f PA.* v) (v PA.+ f PA.* calculateForce params bs b)) $ boids

argsSpec = [
    Arg 0 (Just 's') (Just "seed")
        (argDataOptional "seed" ArgtypeInt) "seed for random number generator",
    Arg 1 (Just 'n') (Just "num")
        (argDataDefaulted "num" ArgtypeInt 10) "number of boids per side",
    Arg 2 Nothing (Just "spacing")
        (argDataDefaulted "spacing" ArgtypeInt 100) "initial spacing of boids",
    Arg 3 Nothing (Just "nearby")
        (argDataDefaulted "radius" ArgtypeFloat 150) "radius that counts as \"nearby\"",
    Arg 4 Nothing (Just "sep")
        (argDataDefaulted "strength" ArgtypeFloat 50000) "strength of separation",
    Arg 5 Nothing (Just "coh")
        (argDataDefaulted "strength" ArgtypeFloat 2) "strength of cohesion",
    Arg 6 Nothing (Just "max_force")
        (argDataDefaulted "force" ArgtypeFloat 100) "maximum force",
    Arg 7 Nothing (Just "max_align")
        (argDataDefaulted "force" ArgtypeFloat 0.1) "maximum alignment force",
    Arg 8 Nothing (Just "width")
        (argDataDefaulted "width" ArgtypeInt 600) "screen width",
    Arg 9 Nothing (Just "height")
        (argDataDefaulted "height" ArgtypeInt 600) "screen Height"
    ]

data SimulationParams = SimulationParams {
    nearbyRadius :: Float, separationStrength :: Float, cohesionStrength :: Float, maxForce :: Float, maxAlignmentSpeed :: Float
}

main :: IO ()
main = do
    args <- parseArgsIO ArgsComplete argsSpec
    seed <- case getArg args 0 of
        Just s -> return s
        Nothing -> randomIO :: IO Int
    let numBoidsPerSide = getRequiredArg args 1
        boidInitialSpacing = getRequiredArg args 2
        screenWidth = getRequiredArg args 8
        screenHeight = getRequiredArg args 9
        params = SimulationParams {
            nearbyRadius = getRequiredArg args 3,
            separationStrength = getRequiredArg args 4,
            cohesionStrength = getRequiredArg args 5,
            maxForce = getRequiredArg args 6, maxAlignmentSpeed = getRequiredArg args 7
            }
    simulate (InWindow "Boids" (screenWidth, screenHeight) (10, 10))
        white 60 (model numBoidsPerSide boidInitialSpacing $ mkStdGen seed) drawModel (stepModel params)
