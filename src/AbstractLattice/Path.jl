"""
     abstract type AbstractPath

Abstract supertype of all paths. A path determines the order to map a lattice to a 1d chain.  
"""
abstract type AbstractPath end

"""
     struct SnakePath <: AbstractPath
"""
struct SnakePath <: AbstractPath end

"""
     struct ZigzagPath <: AbstractPath
"""
struct ZigzagPath <: AbstractPath end

"""
     struct DiagonalPath <: AbstractPath
"""
struct DiagonalPath <: AbstractPath end
