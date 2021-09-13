from django.apps import apps
import os


if __name__ == '__main__':
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'DESIRE.settings')

    for app in apps.get_app_configs():
        print(app.verbose_name, ":")
        for model in app.get_models():
            print("\t", model)