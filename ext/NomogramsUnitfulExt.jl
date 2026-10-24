module NomogramsUnitfulExt

# Units for Nomograms: a variable's `unit` is a Unitful unit, values going in may be
# quantities in any unit of the right dimension, and values coming out carry the
# variable's unit. The core package never sees a quantity -- it works on plain numbers in
# each variable's own unit -- so all of this is the four hooks below.

using Nomograms
using Unitful: Unitful, AbstractQuantity, Units, DimensionError, dimension, unit, ustrip

function Nomograms._strip(u::Units, x::AbstractQuantity, name::Symbol)
    try
        return Float64(ustrip(u, x))
    catch err
        err isa DimensionError || rethrow()
        throw(
            ArgumentError(
                ":$name is in $u, a unit of $(dimension(u)), so $x cannot be " *
                "converted to it",
            ),
        )
    end
end

# A quantity for a variable with no unit has nowhere to go: it would be silently read in
# whatever unit it happened to arrive in.
Nomograms._strip(::Nothing, x::AbstractQuantity, name::Symbol) = throw(
    ArgumentError(
        ":$name has no unit, so the quantity $x has no unit to be converted to. Give the " *
        "variable one -- `unit = u\"$(unit(x))\"`, or a domain in quantities -- or pass a " *
        "plain number.",
    ),
)

Nomograms._attach(u::Units, x::Real) = x * u
Nomograms._infer_unit(x::AbstractQuantity) = unit(x)
Nomograms._unitstring(u::Units) = string(u)

end
