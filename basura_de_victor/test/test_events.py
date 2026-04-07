from src.events import Event, EventsCollection


# --- Event ---

def test_event_creation():
    e = Event(kind="birth", clone_id=(0, 1), rate=0.5)
    assert e.kind == "birth"
    assert e.clone_id == (0, 1)
    assert e.rate == 0.5


def test_event_is_dataclass():
    e1 = Event(kind="death", clone_id=(1,), rate=0.2)
    e2 = Event(kind="death", clone_id=(1,), rate=0.2)
    assert e1 == e2


def test_event_clone_id_is_tuple():
    e = Event(kind="mutation", clone_id=(0, 1, 2), rate=0.1)
    assert isinstance(e.clone_id, tuple)


# --- EventCollection ---

def test_event_collection_starts_empty():
    ec = EventsCollection()
    assert len(ec.events) == 0


def test_add_event():
    ec = EventsCollection()
    e = Event(kind="birth", clone_id=(0,), rate=0.5)
    ec.add_event(e)
    assert len(ec.events) == 1
    assert ec.events[0] == e
    assert ec.get_total_rate() == 0.5


def test_add_multiple_events():
    ec = EventsCollection()
    events = [
        Event(kind="birth", clone_id=(0,), rate=0.5),
        Event(kind="death", clone_id=(1,), rate=0.2),
        Event(kind="mutation", clone_id=(0, 1), rate=0.01),
    ]
    for e in events:
        ec.add_event(e)
    assert len(ec.events) == 3
    assert ec.events == events
    assert ec.get_total_rate() == 0.5 + 0.2 + 0.01


def test_clear_events():
    ec = EventsCollection()
    ec.add_event(Event(kind="birth", clone_id=(0,), rate=0.5))
    ec.add_event(Event(kind="death", clone_id=(1,), rate=0.2))
    ec.clear()
    assert len(ec.events) == 0
    assert ec.get_total_rate() == 0.0


def test_clear_empty_collection():
    ec = EventsCollection()
    ec.clear()  # should not raise
    assert len(ec.events) == 0
    assert ec.get_total_rate() == 0.0