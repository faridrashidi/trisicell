import trisicell as tsc


class TestLogging:
    def test_logging(self):
        tsc.logg.print_version()
        tsc.logg.debug("DEBUG")
        tsc.logg.hint("HINT")
        tsc.logg.info("INFO")
        tsc.logg.warn("WARN")
        assert True
