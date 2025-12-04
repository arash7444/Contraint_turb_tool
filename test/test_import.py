def test_can_import_package():
    import constraint_turb_tool

    assert hasattr(constraint_turb_tool, "generate_mann_fields")
