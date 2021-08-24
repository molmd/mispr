from pymongo import MongoClient


class Connect(object):
    @staticmethod
    def get_connection():
        return MongoClient(
            "mongodb://superuser:idlewide@localhost:27017/?authSource=admin"
        )


if __name__ == "__main__":
    # Connect to local database
    connection = Connect.get_connection()
