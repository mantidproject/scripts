from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer


class AlfFileEventHandler(FileSystemEventHandler):
    def __init__(self):
        self.updated = set()
        path = "C:\Instrument\FileWatcherTestBed"
        observer = Observer()
        observer.schedule(self, path, recursive=True)
        observer.start()
#        try:
#            while True:
#                time.sleep(1)
#        except KeyboardInterrupt:
#            print "hi"
#            observer.stop()
#        observer.join()

    @staticmethod
    def extract_file_name(path):
        split = str(path).split("\\")
        file_name = split[len(split)-1]
        return file_name

    def on_created(self, event):
        file_name = self.extract_file_name(event.src_path)
        if file_name.endswith(".raw"):
            self.updated.add(file_name)
            print "Added: " + file_name

    def on_modified(self, event):
        file_name = self.extract_file_name(event.src_path)
        if file_name.endswith(".raw"):
            self.updated.add(file_name)
            print "Modified: " + file_name

    def on_deleted(self, event):
        file_name = self.extract_file_name(event.src_path)
        if file_name.endswith(".raw"):
            self.updated.discard(file_name)
            print "Removed: " + file_name

    def get_updated(self):
        files = self.updated.copy()
        self.updated.clear()
        return files
