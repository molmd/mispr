from first_connection import Connect

if __name__ == "__main__":
    # Connect to local database
    connection = Connect.get_connection()

    # Access the test database
    db = connection.test
    # Insert a single document
    db.inventory.insert_one(
        {
            "item": "jacket",
            "qty": 48,
            "tags": ["polyester"],
            "size": {"s": 12, "m": 12, "l": 12, "xl": 12},
        }
    )
