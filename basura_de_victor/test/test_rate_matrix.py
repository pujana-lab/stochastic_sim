from src.rate_matrix import RateMatrix
from src.event import Event


# --- EventCollection ---
class TestRateMatrix:
    def test_rate_matrix_starts_empty(self):
        rm = RateMatrix()
        assert len(rm.events) == 0


    def test_add_event(self):
        rm = RateMatrix()
        e = Event(kind="birth", clone_id=(0,), rate=0.5)
        rm.add_event(e)
        assert len(rm.events) == 1
        assert rm.events[0] == e



    def test_add_multiple_events(self):
        rm = RateMatrix()
        events = [
            Event(kind="birth", clone_id=(0,), rate=0.5),
            Event(kind="death", clone_id=(1,), rate=0.2),
            Event(kind="mutation", clone_id=(0, 1), rate=0.01),
        ]
        for e in events:
            rm.add_event(e)
        assert len(rm.events) == 3
        assert rm.events == events
        assert rm.get_total_rate() == 0.5 + 0.2 + 0.01


    def test_clear_events(self):
        rm = RateMatrix()
        rm.add_event(Event(kind="birth", clone_id=(0,), rate=0.5))
        rm.add_event(Event(kind="death", clone_id=(1,), rate=0.2))
        rm.clear()
        assert len(rm.events) == 0
        assert rm.get_total_rate() == 0.0

    def test_clear_empty_collection(self):
        rm = RateMatrix()
        rm.clear()  # should not raise
        assert len(rm.events) == 0
        assert rm.get_total_rate() == 0.0


    def test_get_total_rate_empty(self):
        rm = RateMatrix()
        assert rm.get_total_rate() == 0.0

    def test_get_total_rate_non_empty(self):
        rm = RateMatrix()
        rm.add_event(Event(kind="birth", clone_id=(0,), rate=0.5))
        rm.add_event(Event(kind="death", clone_id=(1,), rate=0.2))
        assert rm.get_total_rate() == 0.7


    def test_choose_event(self):
        """
        Tests choose_event using the inverse-transform method.

        Case A – rates sum to 1.0: [0.5, 0.2, 0.3]
          cumulative boundaries: [0.5, 0.7, 1.0]

        Case B – rates do NOT sum to 1.0: [0.5, 0.25, 0.1]  (total = 0.85)
          normalised:  [0.5882, 0.2941, 0.1176]
          cumulative:  [0.5882, 0.8824, 1.0]
        """
        # ── Case A: normalised rates ──────────────────────────────────────────
        rm = RateMatrix()
        e1 = Event(kind="birth",    clone_id=(0,),   rate=0.5)
        e2 = Event(kind="death",    clone_id=(1,),   rate=0.2)
        e3 = Event(kind="mutation", clone_id=(0, 1), rate=0.3)
        rm.add_event(e1); rm.add_event(e2); rm.add_event(e3)

        assert rm.choose_event(0.0)          == e1   # below first boundary
        assert rm.choose_event(0.3)          == e1   # midpoint of e1 bin
        assert rm.choose_event(0.5 - 1e-9)  == e1   # just before boundary
        assert rm.choose_event(0.5)          == e2   # at boundary → next bin (side='right')
        assert rm.choose_event(0.6)          == e2   # midpoint of e2 bin
        assert rm.choose_event(0.7 - 1e-9)  == e2   # just before second boundary
        assert rm.choose_event(0.7)          == e3   # at boundary → next bin
        assert rm.choose_event(0.99)         == e3   # near 1
        assert rm.choose_event(1.0)          == e3   # edge case: u == 1.0

        # ── Case B: unnormalised rates ────────────────────────────────────────
        rm2 = RateMatrix()
        f1 = Event(kind="birth",    clone_id=(0,),   rate=0.50)
        f2 = Event(kind="death",    clone_id=(1,),   rate=0.25)
        f3 = Event(kind="mutation", clone_id=(0, 1), rate=0.10)
        rm2.add_event(f1); rm2.add_event(f2); rm2.add_event(f3)

        b1 = 0.50 / 0.85   # ≈ 0.5882
        b2 = 0.75 / 0.85   # ≈ 0.8824

        assert rm2.choose_event(0.0)       == f1   # below first boundary
        assert rm2.choose_event(b1 / 2)    == f1   # midpoint of f1 bin
        assert rm2.choose_event(b1 - 1e-9) == f1   # just before boundary
        assert rm2.choose_event(b1)        == f2   # at boundary → next bin
        assert rm2.choose_event(b2 - 1e-9) == f2   # just before second boundary
        assert rm2.choose_event(b2)        == f3   # at boundary → next bin
        assert rm2.choose_event(0.9999)    == f3   # near 1
        assert rm2.choose_event(1.0)       == f3   # edge case: u == 1.0