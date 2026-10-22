@testitem "nicestep rounds up to 1, 2, 2.5 or 5 times a power of ten" tags = [:unit] begin
    @test Nomograms.nicestep(1.0) == 1.0
    @test Nomograms.nicestep(1.4) == 2.0
    @test Nomograms.nicestep(2.3) == 2.5
    @test Nomograms.nicestep(4.0) == 5.0
    @test Nomograms.nicestep(7.0) == 10.0
    @test Nomograms.nicestep(0.013) == 0.02
    @test Nomograms.nicestep(3400.0) == 5000.0
    @test_throws ArgumentError Nomograms.nicestep(0.0)
end

@testitem "linear ticks are round numbers inside the domain" tags = [:unit] begin
    majors, minors = Nomograms.linear_ticks(0.0, 10.0)
    @test sort(majors) == collect(0.0:2.0:10.0)
    @test all(0 .<= minors .<= 10)
    @test isempty(intersect(majors, minors))
    # Subdivisions land exactly, so they print as 0.4 rather than 0.4000000000000001.
    @test 0.4 in minors

    majors, _ = Nomograms.linear_ticks(0.3, 0.7)
    @test all(0.3 .<= majors .<= 0.7)
end

@testitem "linear majors are offered roundest first" tags = [:unit] begin
    majors, _ = Nomograms.linear_ticks(0.0, 10.0)
    # 0 and 10 are multiples of the largest nice step and must be tried before 2 or 6.
    @test findfirst(==(10.0), majors) < findfirst(==(6.0), majors)
    @test findfirst(==(0.0), majors) < findfirst(==(6.0), majors)
end

@testitem "log ticks follow decades" tags = [:unit] begin
    majors, minors = Nomograms.log_ticks(1.0, 1e4)
    @test sort(majors) == [1.0, 10.0, 100.0, 1000.0, 10_000.0]
    @test 20.0 in minors
    @test !(20.0 in majors)

    # Over a short range every mantissa is a major, decades first.
    majors, _ = Nomograms.log_ticks(1.0, 100.0)
    @test 1.0 in majors && 9.0 in majors && 100.0 in majors
    @test findfirst(==(10.0), majors) < findfirst(==(9.0), majors)
    @test findfirst(==(50.0), majors) < findfirst(==(30.0), majors)

    @test_throws ArgumentError Nomograms.log_ticks(0.0, 10.0)
end

@testitem "candidate ticks honour explicit values and order" tags = [:unit] begin
    v = Variable(:a; domain = (0, 10), ticks = [10.0, 0.0, 5.0, 99.0])
    majors, minors = Nomograms.candidate_ticks(v, (0.0, 10.0))
    # Out-of-domain values are dropped; the rest keep the order they were given in.
    @test majors == [10.0, 0.0, 5.0]
    @test isempty(minors)

    vlog = Variable(:a, LogScale(); domain = (1, 1000))
    majors, _ = Nomograms.candidate_ticks(vlog, (1.0, 1000.0))
    @test sort(majors) == [1.0, 10.0, 100.0, 1000.0]
end

@testitem "default_tickformat is compact" tags = [:unit] begin
    @test Nomograms.default_tickformat(10.0) == "10"
    @test Nomograms.default_tickformat(-3.0) == "-3"
    @test Nomograms.default_tickformat(0.5) == "0.5"
    @test Nomograms.default_tickformat(2.4) == "2.4"
    @test Nomograms.default_tickformat(1 / 3) == "0.3333"
end
