from first_connection import Connect

if __name__ == "__main__":
    # Connect to local database
    connection = Connect.get_connection()

    # Access the test database
    db = connection.test

    # Delete single document
    db.inventory.delete_one({"status": "D"})

    # Delete multiple documents
    db.inventory.delete_many({"status": "A"})