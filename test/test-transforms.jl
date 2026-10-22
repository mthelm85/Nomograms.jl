@testitem "transforms round-trip" tags = [:unit] begin
    cases = [
        (LinearScale(), 3.7),
        (LinearScale(2.5, -1.0), 3.7),
        (LogScale(), 8.25),
        (LogScale(2.0), 8.25),
        (ReciprocalScale(), 4.0),
        (PowerScale(2.0), 3.0),
        (PowerScale(-1.5), 3.0),
        (FunctionScale(x -> expm1(x), y -> log1p(y)), 1.3),
        (ComposedScale(LogScale(), PowerScale(2.0)), 3.0),
    ]
    for (t, x) in cases
        y = transform(t, x)
        @test invtransform(t, y) ≈ x
        @test Nomograms.invert(t, y, (x / 10, x * 10)) ≈ x
    end
end

@testitem "transforms without an inverse fall back to root finding" tags = [:unit] begin
    t = FunctionScale(x -> x^3 + x)
    @test !Nomograms.hasinverse(t)
    @test_throws MethodError invtransform(t, 1.0)
    @test Nomograms.invert(t, 10.0, (0.0, 5.0)) ≈ 2.0 atol = 1e-8
    # Outside any bracket, the failure explains itself rather than returning nonsense.
    @test_throws ArgumentError Nomograms.invert(t, 1e6, (0.0, 5.0))
end

@testitem "tick style follows the outermost transform" tags = [:unit] begin
    @test Nomograms.tickstyle(LinearScale()) === :linear
    @test Nomograms.tickstyle(LogScale()) === :log
    @test Nomograms.tickstyle(ReciprocalScale()) === :linear
    @test Nomograms.tickstyle(ComposedScale(LogScale(), ReciprocalScale())) === :log
    @test Nomograms.tickstyle(ComposedScale(LinearScale(), LogScale())) === :linear
end

@testitem "frange reports the span and rejects bad domains" tags = [:unit] begin
    @test Nomograms.frange(LinearScale(), (1.0, 5.0)) == (1.0, 5.0)
    @test Nomograms.frange(ReciprocalScale(), (1.0, 10.0)) == (0.1, 1.0)
    @test all(Nomograms.frange(LogScale(), (1.0, 100.0)) .≈ (0.0, 2.0))

    # A fold-over would map two values to one position: not an axis.
    @test_throws ArgumentError Nomograms.frange(FunctionScale(x -> x^2, sqrt), (-2.0, 2.0))
    # A singularity inside the domain is caught rather than drawn.
    @test_throws ArgumentError Nomograms.frange(ReciprocalScale(), (-1.0, 1.0))
    @test_throws ArgumentError Nomograms.frange(LinearScale(), (5.0, 1.0))
end
