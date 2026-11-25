# Number lists of the text fields (Python: tests/test_gui.py test_number_lists,
# test_number_list_syntax, test_number_lists_in_display_units,
# test_points_and_formatting; tests/test_case_tab.py parse_points in units).
using HTGEquilibrium, Test

@testset "number lists" begin
    @testset "lists and ranges" begin
        @test parse_numbers("300, 350;400 450") == [300, 350, 400, 450]
        @test parse_numbers("300:400:25") == [300, 325, 350, 375, 400]
        @test parse_numbers("300:390:25") == [300, 325, 350, 375]
        @test parse_numbers("100, 300:350:50") == [100, 300, 350]
        @test parse_numbers("") == [] && parse_numbers(nothing) == []
        @test_throws NumberListError parse_numbers("300:400:0")
        @test_throws NumberListError parse_numbers("hot")
    end

    @testset "list syntax (spec WP1 item 4)" begin
        @test parse_numbers("300:400 n=5") == [300, 325, 350, 375, 400]
        @test parse_numbers("300, 350, 400") == [300, 350, 400]
        @test parse_numbers("300:400:50") == [300, 350, 400]
        @test parse_numbers("250, 300:400:50") == [250, 300, 350, 400]
        @test parse_numbers("400; 250\t300:400 n=3") == [250, 300, 350, 400]   # sorted, duplicates removed
        @test parse_numbers("300 : 400  n = 3") == [300, 350, 400]
        @test parse_numbers("0:1:0.1") == [i / 10 for i in 0:10]
        @test parse_numbers("350 300 300"; sort=false) == [350, 300, 300]
        for (text, message) in [("n=5", "must follow a range"), ("300:400 n=1", "2 or more"),
                                ("300:400 n=2.5", "whole number"), ("300:400", "needs a step or a count"),
                                ("300:400:50 n=3", "both a step and a count"), ("400:300:50", "start is above"),
                                ("300:400:-5", "above zero"), ("1:2:3:4", "not a range"), ("nan", "not a number"),
                                ("0:1e9:1", "more than"), ("300:400 n=20000", "at most")]
            err = try
                parse_numbers(text)
                nothing
            catch e
                e
            end
            @test err isa NumberListError && err isa HTGEquilibrium.ValueError && occursin(message, err.msg)
        end
        @test sprint(showerror, NumberListError("x")) == "x"
        @test (try parse_numbers("n=5") catch e e.msg end) == "\"n=5\" must follow a range start:stop, e.g. 300:400 n=5."
        @test (try parse_numbers("300:400 n=1") catch e e.msg end) ==
              "\"n=1\": the count after n= must be a whole number of 2 or more."
        @test (try parse_numbers("hot") catch e e.msg end) ==
              "\"hot\" is not a number. Write numbers like 300, 350, a range 300:400:50 or 300:400 n=5."
    end

    @testset "display units" begin
        @test parse_numbers("573.15:673.15 n=3"; unit="K") == [300, 350, 400]
        @test parse_numbers("24, 30"; unit="MPa") == [240, 300]
        @test parse_numbers("2.4e7"; unit="Pa") == [240]
        @test parse_numbers("1"; unit="atm") == [1.01325]
        @test parse_numbers("300"; unit="C") == [300] && parse_numbers("300"; unit="°C") == [300]
        err = try parse_numbers("300"; unit="F") catch e e end
        @test err isa Units.UnitError && startswith(err.msg, "Unknown unit")
    end

    @testset "points and formatting" begin
        @test parse_points("300, 240\n\n350 300\n") == [[300, 240], [350, 300]]
        @test_throws NumberListError parse_points("300")
        @test (try parse_points("300") catch e e.msg end) == "\"300\": give a temperature and a pressure."
        @test parse_points("573.15, 24\n623.15 30", "K", "MPa") == [[300, 240], [350, 300]]
        @test parse_points("") == []
        @test format_numbers([300, 350.5]) == "300, 350.5"
        @test format_numbers(240, "Pa") == "2.4e7"
        @test format_numbers([300], "K") == "573.15"
        @test format_numbers(nothing) == ""
        @test HTGEquilibrium.clean_number(299.99999999999994) == 300
    end
end
