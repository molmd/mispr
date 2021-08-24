from first_connection import Connect

if __name__ == "__main__":
    # Connect to local database
    connection = Connect.get_connection()

    # Access the test database
    db = connection.test

    # Update a single document in the inventory collection
    db.inventory.update_one(
        {"item": "paper"},
        {
            "$set": {"size.uom": "cm", "status": "P"},
            "$currentDate": {"lastModified": True},
        },
    )

    # Update multiple documents
    db.inventory.update_many(
        {"qty": {"$lt": 50}},
        {
            "$set": {"size.uom": "in", "status": "P"},
            "$currentDate": {"lastModified": True},
        },
    )
