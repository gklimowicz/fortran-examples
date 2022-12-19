import whizard

# Imitate API_UT_CC_1

wz = whizard.Whizard()
wz.option("logfile", "whizard_1_py.log")
wz.init()
del(wz)

# API_UT_CC_2
wz = whizard.Whizard()
wz.option("logfile", "whizard_2_py.log")
wz.option("job_id", "whizard_2_py_ID")
wz.init()

try:
    job_id = wz.get_string("$job_id")
except IndexError as e:
    print(str(e))
else:
    print(f"$job_id = {job_id}")

wz.set_double("sqrts", 100)
wz.set_int("n_events", 3)
wz.set_bool("?unweighted", False)
wz.set_string("$sample", "foobar")

assert(wz.get_double("sqrts") == 100.)
assert(wz.get_int("n_events") == 3)
assert(not wz.get_bool("?unweighted"))
assert(wz.get_string("$sample") == "foobar")

wz.set_bool("?unweighted", True)
assert(wz.get_bool("?unweighted"))
del(wz)

# API_UT_CC_3
wz = whizard.Whizard()
wz.option("logfile", "whizard_1_py_3.log")
wz.option("model", "QED")
wz.init()

assert(wz.flv_string(11) == '"e-"')
assert(wz.flv_array_string([11, 13, 15]) == '"e-":"m-":"t-"')
del(wz)

# API UT CC 4
wz = whizard.Whizard()
wz.option("logfile", "whizard_1_py_4.log")
wz.option("library", "whizard_1_py_4_lib")
wz.option("model", "QED")
wz.option("rebuild", "true")
wz.init()

wz.command("process whizard_1_py_4_p = e1, E1 => e2, E2")
wz.command("sqrts = 10")
wz.command("iterations = 1:100")
wz.set_int("seed", 0)
integral, error = wz.get_integration_result("whizard_1_py_4_p")
print(integral, error)

wz.command("integrate (whizard_1_py_4_p)")
integral, error = wz.get_integration_result("whizard_1_py_4_p")
sqrts = wz.get_double("sqrts")
print(f"sqrts= {sqrts:5.1f} GeV")
print(f"sigma ={integral:5.1f} pb")
print(f"error= {error:5.1f} pb")


